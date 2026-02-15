/**
 * @file   grt_modsum.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 *     根据计算的相速度频散(本征值)，使用模态叠加（modal summation）计算面波格林函数
 * 
 */

#include "grt.h"

#include "grt/common/matrix.h"
#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/util.h"
#include "grt/common/search.h"
#include "grt/common/myfftw.h"
#include "grt/common/sacio2.h"

#include "grt/dynamic/source.h"
#include "grt/dynamic/layer.h"
#include "grt/common/travt.h"
#include "grt/modal/modal_util.h"
#include "grt/modal/secular.h"
#include "grt/modal/eigenfn.h"
#include "grt/modal/energy.h"


// 一些变量的非零默认值
#define GRT_GREENFN_N_UPSAMPLE    1
#define GRT_GREENFN_G_EX       true
#define GRT_GREENFN_G_VF       true
#define GRT_GREENFN_G_HF       true
#define GRT_GREENFN_G_DC       true

/** 该子模块的参数控制结构体 */
typedef struct {
    /* 输入频散文件 */
    struct {
        bool active;
        char *s_filepath;
    } C;
    /* 阶数范围 */
    struct {
        bool active;
        size_t *modes;
        size_t nmode;
    } N;
    /* 频率范围 */
    struct {
        bool active;
        size_t nf;
        real_t *freqs;
    } F;
    /* 输入震源和台站深度，计算面波格林函数 */
    struct {
        bool active;
        real_t depsrc;
        real_t deprcv;
        char *s_depsrc;
        char *s_deprcv;
    } D;
    /* 面波格林函数输出路径 */
    struct {
        bool active;
        char *s_output_dir;
    } O;
    /** 时间延迟 */
    struct {
        bool active;
        real_t delayT0;
        real_t delayV0;
        bool refFirstP;  ///< 是否参考初至P
    } E;
    /* 震中距数组 */
    struct {
        bool active;
        char *s_raw;
        char **s_rs;
        real_t *rs;
        size_t nr;
    } R;
    /** 多线程 */
    struct {
        bool active;
        int nthreads; ///< 线程数
    } P;
    /** 输出哪些震源的格林函数 */
    struct {
        bool active;
        bool doEX;
        bool doVF;
        bool doHF;
        bool doDC;
    } G;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;
    /* 升采样倍数 */
    struct {
        bool active;
        int upsample_n;
    } W;

} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // C
    GRT_SAFE_FREE_PTR(Ctrl->C.s_filepath);
    // N
    GRT_SAFE_FREE_PTR(Ctrl->N.modes);

    // D
    GRT_SAFE_FREE_PTR(Ctrl->D.s_depsrc);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_deprcv);
    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
    // R
    GRT_SAFE_FREE_PTR_ARRAY(Ctrl->R.s_rs, Ctrl->R.nr);
    GRT_SAFE_FREE_PTR(Ctrl->R.s_raw);
    GRT_SAFE_FREE_PTR(Ctrl->R.rs);

    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
    printf("[WAITING TO FINISH]\n");
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 先为个别参数设置非0初始值
    Ctrl->W.upsample_n = GRT_GREENFN_N_UPSAMPLE;

    Ctrl->G.doEX = GRT_GREENFN_G_EX;
    Ctrl->G.doVF = GRT_GREENFN_G_VF;
    Ctrl->G.doHF = GRT_GREENFN_G_HF;
    Ctrl->G.doDC = GRT_GREENFN_G_DC;

    int opt;
    while ((opt = getopt(argc, argv, ":C:N::F:D:O:E:R:G:W:P:esh" )) != -1) {
        switch (opt) {
            // 频散文件
            case 'C':
                Ctrl->C.active = true;
                Ctrl->C.s_filepath = strdup(optarg);
                break;

            // 阶数范围 -N[n1][/n2][/dn]
            case 'N':
                Ctrl->N.active = true;
                // 没有参数，直接跳过
                if(optarg==NULL)  break;
                {
                    ssize_t a1, a2, dif;
                    a1 = a2 = 0;
                    dif = 1;
                    int nscan = sscanf(optarg, "%zd/%zd/%zd", &a1, &a2, &dif);
                    if( !(nscan == 1 || nscan == 2 || nscan == 3)){
                        GRTBadOptionError(N, "");
                    }
                    if(nscan == 1 && a1 < 0){
                        GRTBadOptionError(N, "Can't set a single negative value.");
                    }
                    if(nscan == 1){
                        a2 = a1;
                        dif = 1;
                    }
                    else if(nscan >=2){
                        if(dif <= 0){
                            GRTBadOptionError(N, "Can't set negative dn(%zd).", dif);
                        }
                        if(a1 < 0 || a2 < 0){
                            GRTBadOptionError(N, "Can't set nonpositive n1(%zd), n2(%zd).", a1, a2);
                        }
                        if(a1 > a2){
                            GRTBadOptionError(N, "n1(%zd) > n2(%zd).", a1, a2);
                        }
                    }
                    Ctrl->N.nmode = floor((a2-a1)/dif) + 1;
                    Ctrl->N.modes = (size_t *)calloc(Ctrl->N.nmode, sizeof(size_t));
                    for(size_t i=0; i < Ctrl->N.nmode; ++i){
                        Ctrl->N.modes[i] = a1 + dif*i;
                    }
                }
                break;

            // 频率值 -Ff1/f2  指定频率范围
            case 'F':
                Ctrl->F.active = true;
                Ctrl->F.freqs = (real_t *)calloc(2, sizeof(real_t));
                if(2 != sscanf(optarg, "%lf/%lf", Ctrl->F.freqs, Ctrl->F.freqs+1)){
                    GRTBadOptionError(F, "");
                };
                if(Ctrl->F.freqs[0] < 0.0){
                    GRTBadOptionError(F, "Can't set negative freq1(%f).", Ctrl->F.freqs[0]);
                }
                if(Ctrl->F.freqs[1] < 0.0){
                    GRTBadOptionError(F, "Can't set negative freq2(%f).", Ctrl->F.freqs[1]);
                }
                if(Ctrl->F.freqs[0] > Ctrl->F.freqs[1]){
                    GRTBadOptionError(F, "Positive freq1 should be less than positive freq2.");
                }
                break;

            // 输入震源和台站深度、震中距和保存目录，计算面波格林函数
            // 震源和场点深度， -Ddepsrc/deprcv
            case 'D':
                Ctrl->D.active = true;
                Ctrl->D.s_depsrc = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                Ctrl->D.s_deprcv = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                if(2 != sscanf(optarg, "%[^/]/%s", Ctrl->D.s_depsrc, Ctrl->D.s_deprcv)){
                    GRTBadOptionError(D, "");
                };
                if(1 != sscanf(Ctrl->D.s_depsrc, "%lf", &Ctrl->D.depsrc)){
                    GRTBadOptionError(D, "");
                }
                if(1 != sscanf(Ctrl->D.s_deprcv, "%lf", &Ctrl->D.deprcv)){
                    GRTBadOptionError(D, "");
                }
                if(Ctrl->D.depsrc < 0.0 || Ctrl->D.deprcv < 0.0){
                    GRTBadOptionError(D, "Negative value in -D is not supported.");
                }
                break;
            
            // 输出路径 -Ooutput_dir
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_output_dir = strdup(optarg);
                break;

            // 时间延迟 -E[p]<t0>[/<v0>]
            case 'E':
                Ctrl->E.active = true;
                if(optarg[0] == 'p'){
                    if(1 != sscanf(optarg+1, "%lf", &Ctrl->E.delayT0)){
                        GRTBadOptionError(E, "");
                    };
                    if(Ctrl->E.delayT0 >= 0.0){
                        GRTBadOptionError(E, "Can't set positive t0(%f) in -Ep.", Ctrl->E.delayT0);
                    }
                    Ctrl->E.refFirstP = true;
                } else {
                    if(0 == sscanf(optarg+1, "%lf/%lf", &Ctrl->E.delayT0, &Ctrl->E.delayV0)){
                        GRTBadOptionError(E, "");
                    };
                    if(Ctrl->E.delayV0 < 0.0){
                        GRTBadOptionError(E, "Can't set negative v0(%f) in -E.", Ctrl->E.delayV0);
                    }
                }
                break;

            // 震中距数组，-R<r1>,<r2>[,...]|<r1>/<r2>/<dr>|<file>
            case 'R':
                Ctrl->R.active = true;
                Ctrl->R.s_raw = strdup(optarg);
                {
                    real_t a1, a2, delta;
                    
                    // 如果输入仅由数字、小数点和间隔符组成，则直接读取
                    if(grt_string_composed_of(optarg, GRT_NUM_STR "eE+-" ".,")){
                        Ctrl->R.s_rs = grt_string_split(optarg, ",", &Ctrl->R.nr);
                    }
                    // 尝试按照 <r1>/<r2>/<dr> 读取
                    else if(3 == sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        if(delta <= 0){
                            GRTBadOptionError(R, "Can't set nonpositive dr(%f)", delta);
                        }
                        if(a1 > a2){
                            GRTBadOptionError(R, "r1(%f) > r2(%f).", a1, a2);
                        }

                        // 根据最大的小数位数来设置 s_rs
                        int place = 0;
                        real_t test[3] = {0};
                        test[0] = a1;
                        test[1] = a2;
                        test[2] = delta;
                        do {
                            bool agree = true;
                            for(size_t i = 0; i < 3; ++i)  {
                                real_t new_frac = test[i] - (int)test[i];
                                test[i] *= 10.0;
                                // 如果乘以10后变成整数（或非常接近整数）
                                if (fabs(new_frac) > 1e-8) {
                                    agree = false;
                                }
                            }
                            if(agree) break;
                            place++;
                        } while (place < 8);  // 最多小数点后 8 位

                        Ctrl->R.nr = floor((a2-a1)/delta) + 1;
                        Ctrl->R.s_rs = (char **)calloc(Ctrl->R.nr, sizeof(char*) * Ctrl->R.nr);
                        for(size_t ir = 0; ir < Ctrl->R.nr; ++ir){
                            GRT_SAFE_ASPRINTF(&Ctrl->R.s_rs[ir], "%.*f", place, a1 + delta*ir);
                        }
                    }
                    // 否则从文件读取
                    else {
                        FILE *fp = GRTCheckOpenFile(optarg, "r");
                        Ctrl->R.s_rs = grt_string_from_file(fp, &Ctrl->R.nr);
                        fclose(fp);
                    }
                        
                    // 转为浮点数
                    Ctrl->R.rs = (real_t*)realloc(Ctrl->R.rs, sizeof(real_t)*(Ctrl->R.nr));
                    for(size_t i=0; i<Ctrl->R.nr; ++i){
                        Ctrl->R.rs[i] = atof(Ctrl->R.s_rs[i]);
                        if(Ctrl->R.rs[i] < 0.0){
                            GRTBadOptionError(R, "Can't set negative epicentral distance(%f).", Ctrl->R.rs[i]);
                        }
                    }
                }
                break;

            // 选择要计算的格林函数 -Ge|v|h|s
            case 'G': 
                Ctrl->G.active = true;
                // 先全部置否
                Ctrl->G.doEX = Ctrl->G.doVF = Ctrl->G.doHF = Ctrl->G.doDC = false;
                for(size_t i=0; i < strlen(optarg); ++i){
                    switch (optarg[i]) {
                        case 'e':
                            Ctrl->G.doEX = true;
                            break;
                        case 'v':
                            Ctrl->G.doVF = true;
                            break;
                        case 'h':
                            Ctrl->G.doHF = true;
                            break;
                        case 's':
                            Ctrl->G.doDC = true;
                            break;
                        default:
                            GRTBadOptionError(G, "unknown type %c.", optarg[i]);
                            break;
                    }
                }
                // 至少要有一个真
                if(!(Ctrl->G.doEX || Ctrl->G.doVF || Ctrl->G.doHF || Ctrl->G.doDC)){
                    GRTBadOptionError(G, "At least set one true value.");
                }
                break;

            // 升采样倍数
            case 'W':
                Ctrl->W.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->W.upsample_n)){
                    GRTBadOptionError(W, "");
                };
                if(Ctrl->W.upsample_n <= 0){
                    GRTBadOptionError(W, "need positive integer, but get (%d).", Ctrl->W.upsample_n);
                }
                break;
            
            // 多线程数 -Pnthreads
            case 'P':
                Ctrl->P.active = true;
                if(1 != sscanf(optarg, "%d", &Ctrl->P.nthreads)){
                    GRTBadOptionError(P, "");
                };
                if(Ctrl->P.nthreads <= 0){
                    GRTBadOptionError(P, "Nonpositive value is not supported.");
                }
                grt_set_num_threads(Ctrl->P.nthreads);
                break;

            // 是否计算位移空间导数
            case 'e':
                Ctrl->e.active = true;
                break;

            GRT_Common_Options_in_Switch((char)(optopt));

        }
    }


    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, C);
    GRTCheckOptionActive(Ctrl, D);
    GRTCheckOptionActive(Ctrl, O);
    GRTCheckOptionActive(Ctrl, R);

    // -N 的默认项，仅处理基阶
    if( ! Ctrl->N.active ){
        Ctrl->N.nmode = 1;
        Ctrl->N.modes = (size_t *)calloc(Ctrl->N.nmode, sizeof(size_t));
        Ctrl->N.modes[0] = 0;
    }
    
    // 建立保存目录
    GRTCheckMakeDir(Ctrl->O.s_output_dir);

    // 在目录中保留命令
    char *dummy = NULL;
    GRT_SAFE_ASPRINTF(&dummy, "%s/command", Ctrl->O.s_output_dir);
    FILE *fp = GRTCheckOpenFile(dummy, "a");
    fprintf(fp, GRT_MAIN_COMMAND " ");  // 主程序名
    for(int i=0; i<argc; ++i){
        fprintf(fp, "%s ", argv[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    GRT_SAFE_FREE_PTR(dummy);

}


static void grt_modsum_grn_Rayl(
    const GRT_MODEL1D *mod1d, const cplx_t (*mod_potRaylLove_Down)[GRT_RAYL_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_RAYL_DIM],
    const cplx_t *egyint, const size_t nr, const real_t *rs, const size_t iw, 
    bool calc_upar, pcplxChnlGrid *grn, pcplxChnlGrid *grn_uiz, pcplxChnlGrid *grn_uir)
{
    // 注意模型中并未插入震源层和接收层
    size_t isrc = mod1d->isrc;
    size_t ircv = mod1d->ircv;

    real_t eigenK = mod1d->k;

    // 震源层物性参数
    cplx_t src_mu = mod1d->mu[isrc];
    real_t src_rho = mod1d->Rho[isrc];
    cplx_t src_a = eigenK * mod1d->xa[isrc];
    cplx_t src_b = eigenK * mod1d->xb[isrc];
    cplx_t src_Omg = eigenK*eigenK * (1.0 - 0.5*mod1d->cbcb[isrc]);

    // 计算震源系数
    cplxChnlGrid src_coefD = {0};
    cplxChnlGrid src_coefU = {0};
    grt_source_coef_PSV(mod1d, src_coefD, src_coefU);

    // 构造所需的震源系数
    cplx_t fcoef[GRT_SRC_M_NUM][GRT_RAYL_DIM];
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        fcoef[i][0] =   eigenK*(src_coefD[i][0] - src_coefU[i][0]) - src_b*(src_coefD[i][1] + src_coefU[i][1]);
        fcoef[i][1] = - src_a *(src_coefD[i][0] + src_coefU[i][0]) + eigenK*(src_coefD[i][1] - src_coefU[i][1]);
        fcoef[i][2] =  2.0*src_mu*(src_Omg*(src_coefD[i][0] - src_coefU[i][0]) - eigenK*src_b*(src_coefD[i][1] + src_coefU[i][1]));
        fcoef[i][3] =  2.0*src_mu*(- eigenK*src_a*(src_coefD[i][0] + src_coefU[i][0]) + src_Omg*(src_coefD[i][1] - src_coefU[i][1]));
    }

    // 计算震源处的本征函数
    cplx_t T0[GRT_RAYL_DIM][GRT_RAYL_DIM] = {0};
    cplx_t src_eigenfn[4] = {0};
    grt_get_eigenfn_single_depth_Rayl(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, T0, false, mod1d->depsrc, isrc, src_eigenfn);

    // 构造公式中的 Fm^R 项
    cplx_t F_item[GRT_SRC_M_NUM] = {0};
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        F_item[i] =   fcoef[i][0]*src_eigenfn[3] + fcoef[i][1]*src_eigenfn[2]
                    - fcoef[i][2]*src_eigenfn[1] - fcoef[i][3]*src_eigenfn[0];
    }

    // 计算接收处的本征函数
    cplx_t rcv_eigenfn[4] = {0};
    cplx_t rcv_eigenfn_z[4] = {0};
    grt_get_eigenfn_single_depth_Rayl(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, T0, false, mod1d->deprcv, ircv, rcv_eigenfn);
    if(calc_upar){
        cplx_t xa = mod1d->xa[ircv];
        cplx_t xb = mod1d->xb[ircv];
        bool   isLiquid = mod1d->isLiquid[ircv];
        real_t k = mod1d->k;

        cplx_t ak = k*k*xa;
        cplx_t bk = k*k*xb;
        cplx_t bb = xb*bk;
        cplx_t aa = xa*ak;
        cplx_t T0_z[GRT_RAYL_DIM][GRT_RAYL_DIM] = {
            {ak, bb, -ak, bb},
            {aa, bk, aa, -bk}
        };

        if(! isLiquid){
            grt_get_eigenfn_single_depth_Rayl(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, T0_z, true, mod1d->deprcv, ircv, rcv_eigenfn_z);
        }
    }

    // 计算该频点、该本征值下的频谱
    cplx_t i_m[GRT_MORDER_MAX+1] = {1.0, I, -1.0};
    cplx_t dom = 4.0*egyint[1] - 2.0*egyint[2]/eigenK;
    cplx_t wcoef, qcoef;
    cplx_t wcoefz, qcoefz;
    for(int im = 0; im < GRT_SRC_M_NUM; ++im){
        int modr = GRT_SRC_M_ORDERS[im];
        // Rayleigh 波仅考虑 q, w 项
        
        cplx_t aa = - 1.0/(4.0*src_rho*GRT_SQUARE(mod1d->omega));

        // Z  
        wcoef = rcv_eigenfn[1] * i_m[modr] * F_item[im] / dom * aa;
        wcoefz = rcv_eigenfn_z[1] * i_m[modr] * F_item[im] / dom * aa;
        // R
        qcoef = I * rcv_eigenfn[0] * i_m[modr] * F_item[im] / dom * aa;
        qcoefz = I * rcv_eigenfn_z[0] * i_m[modr] * F_item[im] / dom * aa;

        for(size_t ir = 0; ir < nr; ++ir){
            // 草稿推导有误？似乎应该使用 H^(2)(kr)
            cplx_t ekr = sqrt(2.0 / (PI*eigenK*rs[ir])) * exp( - I * (eigenK*rs[ir] + QUARTERPI));
            grn[ir][im][0][iw] += wcoef * ekr;
            grn[ir][im][1][iw] += - qcoef * ekr;
            if(calc_upar){
                grn_uiz[ir][im][0][iw] += wcoefz * ekr;
                grn_uiz[ir][im][1][iw] += - qcoefz * ekr;
            }
        }
    }

}



static void grt_modsum_grn_Love(
    const GRT_MODEL1D *mod1d, const cplx_t (*mod_potRaylLove_Down)[GRT_LOVE_DIM], const cplx_t (*mod_potRaylLove_Up)[GRT_LOVE_DIM],
    const cplx_t *egyint, const size_t nr, const real_t *rs, const size_t iw, 
    bool calc_upar, pcplxChnlGrid *grn, pcplxChnlGrid *grn_uiz, pcplxChnlGrid *grn_uir)
{
    // 注意模型中并未插入震源层和接收层
    size_t isrc = mod1d->isrc;
    size_t ircv = mod1d->ircv;

    real_t eigenK = mod1d->k;

    // 震源层物性参数
    cplx_t src_mu = mod1d->mu[isrc];
    real_t src_rho = mod1d->Rho[isrc];
    cplx_t src_b = eigenK * mod1d->xb[isrc];

    // 计算震源系数
    cplxChnlGrid src_coefD = {0};
    cplxChnlGrid src_coefU = {0};
    grt_source_coef_SH(mod1d, src_coefD, src_coefU);

    // 构造所需的震源系数
    cplx_t fcoef[GRT_SRC_M_NUM][GRT_LOVE_DIM] = {0};
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        fcoef[i][0] =   eigenK*(src_coefD[i][2] - src_coefU[i][2]);
        fcoef[i][1] = - src_mu*eigenK*src_b*(src_coefD[i][2] + src_coefU[i][2]);
    }

    // 计算震源处的本征函数
    cplx_t T0[GRT_LOVE_DIM][GRT_LOVE_DIM] = {0};
    cplx_t src_eigenfn[4] = {0};
    grt_get_eigenfn_single_depth_Love(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, T0, false, mod1d->depsrc, isrc, src_eigenfn);

    // 构造公式中的 Fm^L 项
    cplx_t F_item[GRT_SRC_M_NUM] = {0};
    for(int i = 0; i < GRT_SRC_M_NUM; ++i){
        F_item[i] = fcoef[i][0]*src_eigenfn[1] - fcoef[i][1]*src_eigenfn[0];
    }

    // 计算接收处的本征函数
    cplx_t rcv_eigenfn[4] = {0};
    cplx_t rcv_eigenfn_z[4] = {0};
    grt_get_eigenfn_single_depth_Love(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, T0, false, mod1d->deprcv, ircv, rcv_eigenfn);
    if(calc_upar){
        cplx_t xb = mod1d->xb[ircv];
        bool   isLiquid = mod1d->isLiquid[ircv];
        real_t k = mod1d->k;

        cplx_t bk = k*k*xb;
        cplx_t T0_z[GRT_LOVE_DIM][GRT_LOVE_DIM] = {{bk, -bk}};

        if(! isLiquid){
            grt_get_eigenfn_single_depth_Love(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, T0_z, true, mod1d->deprcv, ircv, rcv_eigenfn_z);
        }
    }

    // 计算该频点、该本征值下的频谱
    cplx_t i_m[GRT_MORDER_MAX+1] = {1.0, I, -1.0};
    cplx_t dom = 4.0*egyint[1];
    cplx_t vcoef;
    cplx_t vcoefz;
    for(int im = 0; im < GRT_SRC_M_NUM; ++im){
        int modr = GRT_SRC_M_ORDERS[im];
        // Love 波仅考虑 v 项
        if(modr == 0)  continue;
        
        cplx_t aa = - 1.0/(4.0*src_rho*GRT_SQUARE(mod1d->omega));

        // T
        vcoef = I * rcv_eigenfn[0] * i_m[modr] * F_item[im] / dom * aa;
        vcoefz = I * rcv_eigenfn_z[0] * i_m[modr] * F_item[im] / dom * aa;

        for(size_t ir = 0; ir < nr; ++ir){
            // 草稿推导有误？似乎应该使用 H^(2)(kr)
            cplx_t ekr = sqrt(2.0 / (PI*eigenK*rs[ir])) * exp( - I * (eigenK*rs[ir] + QUARTERPI));
            grn[ir][im][2][iw] += vcoef * ekr;
            if(calc_upar){
                grn_uiz[ir][im][2][iw] += vcoefz * ekr;
            }
        }
    }

}

static void grt_modsum_grn(
    const GRT_MODEL1D *mod1d, const DISPER_TYPE wtype, const size_t ncols, 
    const cplx_t (*mod_potRaylLove_Down)[ncols], const cplx_t (*mod_potRaylLove_Up)[ncols],
    const cplx_t *egyint, const size_t nr, const real_t *rs, const size_t iw, 
    bool calc_upar, pcplxChnlGrid *grn, pcplxChnlGrid *grn_uiz, pcplxChnlGrid *grn_uir)
{
    if(wtype == GRT_DISPERSION_RAYL && ncols == GRT_RAYL_DIM){
        grt_modsum_grn_Rayl(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, egyint, nr, rs, iw, calc_upar, grn, grn_uiz, grn_uir);
    } 
    else if(wtype == GRT_DISPERSION_LOVE && ncols == GRT_LOVE_DIM) {
        grt_modsum_grn_Love(mod1d, mod_potRaylLove_Down, mod_potRaylLove_Up, egyint, nr, rs, iw, calc_upar, grn, grn_uiz, grn_uir);
    }
    else {
        GRTRaiseError("Wrong execution.");
    }

    // 对一些特殊情况的修正
    // 当震源和场点均位于地表时，可理论验证DS分量恒为0，这里直接赋0以避免后续的精度干扰
    if(mod1d->depsrc == 0.0 && mod1d->deprcv == 0.0)
    {
        for(size_t ir = 0; ir < nr; ++ir){
            for(int c = 0; c < GRT_CHANNEL_NUM; ++c){
                grn[ir][GRT_SRC_M_DS_INDEX][c][iw] = 0.0;
            }
        }
    }
}



/** 将一条数据反变换回时间域再进行处理，保存到SAC文件 */
static void write_one_to_sac(
    const char *srcname, const char ch, GRT_FFTW_HOLDER *fh, const real_t wI, 
    SACTRACE *sac, const char *s_output_subdir, const char *s_prefix,
    const int sgn, bool skipImagComps, const cplx_t *grncplx)
{
    snprintf(sac->hd.kcmpnm, sizeof(sac->hd.kcmpnm), "%s%s%c", s_prefix, srcname, ch);
    
    char *s_outpath = NULL;
    GRT_SAFE_ASPRINTF(&s_outpath, "%s/%s.sac", s_output_subdir, sac->hd.kcmpnm);

    // 执行fft任务会修改数组，需重新置零
    grt_reset_fftw_holder_zero(fh);
    
    // 赋值复数，包括时移
    cplx_t cfac, ccoef;
    cfac = exp(I*PI2*fh->df*sac->hd.b);
    ccoef = sgn;
    // 只赋值有效长度，其余均为0
    for(size_t i = 0; i < fh->nf_valid; ++i){
        fh->W_f[i] = grncplx[i] * ccoef;
        ccoef *= cfac;
    }


    if(! fh->naive_inv){
        // 发起fft任务 
        fftw_execute(fh->plan);
    } else {
        grt_naive_inverse_transform_double(fh);
    }

    // 归一化，并处理虚频
    // 并转为 SAC 需要的单精度类型
    real_t fac, coef;
    coef = fh->df;
    fac = 1.0;
    if (! skipImagComps) {
        coef *= exp(sac->hd.b * wI);
        fac = exp(wI*fh->dt);
    }
    for(size_t i = 0; i < fh->nt; ++i){
        sac->data[i] = fh->w_t[i] * coef;
        coef *= fac;
    }
    

    // 以sac文件保存到本地
    grt_write_SACTRACE(s_outpath, sac);

    GRT_SAFE_FREE_PTR(s_outpath);
}



/* 子模块主函数 */
int modsum_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 读取频散
    char *modelpath = NULL;
    EIGENV_INFO *eigmet = (EIGENV_INFO *)calloc(1, sizeof(EIGENV_INFO));
    grt_read_cdisp(Ctrl->C.s_filepath, eigmet, &modelpath);

    // 读取模型（不插入震源和台站的虚拟层）
    GRT_MODEL1D *mod1d = NULL;
    if((mod1d = grt_read_mod1d_from_file(modelpath, -1.0, -1.0, true)) ==NULL){
        exit(EXIT_FAILURE);
    }

    // 当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    // 程序结束前会输出对应警告
    if(mod1d->isLiquid[mod1d->isrc]){
        Ctrl->G.doHF = Ctrl->G.doVF = Ctrl->G.doDC = false;
    }

    // 不再插入虚拟层，直接记录震源和台站所在的真实层位
    mod1d->depsrc = Ctrl->D.depsrc;
    mod1d->deprcv = Ctrl->D.deprcv;
    mod1d->ircvup = (Ctrl->D.depsrc >= Ctrl->D.deprcv);
    for(mod1d->isrc=0; mod1d->isrc < mod1d->n-1; ++mod1d->isrc){
        if(mod1d->Dep[mod1d->isrc+1] >= Ctrl->D.depsrc)  break; 
    }
    for(mod1d->ircv=0; mod1d->ircv < mod1d->n-1; ++mod1d->ircv){
        if(mod1d->Dep[mod1d->ircv+1] >= Ctrl->D.deprcv)  break; 
    }

    // 根据命令行参数确定出所需的部分频散信息
    EIGENFN_INFO *eigfnmet = (EIGENFN_INFO *)calloc(1, sizeof(EIGENFN_INFO));
    grt_filter_eigenfn_info(
        Ctrl->F.nf, Ctrl->F.freqs, true,
        Ctrl->N.nmode, Ctrl->N.modes, eigmet, eigfnmet);
    // 使用模型层计算能量积分
    eigfnmet->cpar_nz = mod1d->n;
    
    // 如果要计算面波格林函数，对生成频散文件的参数有限制，在 eigenv 必须以 -F0/fmax/df 的形式使用
    // 不能只传入一个频率值
    if(eigfnmet->nf <= 1){
        GRTRaiseError("Modal summation need at least two sampled freqs, but only found %zu.\n", eigfnmet->nf);
    }
    // 原频散文件中，起始频率必须等于频率间隔，且不能以周期的形式给定
    if(fabs(2.0*eigmet->freqs[0] - eigmet->freqs[1]) > 1e-8){
        GRTRaiseError("Modal summation need dispersion file generated like `grt eigenv -F0/fmax/df`.\n");
    }

    // 初始化 eigfn
    eigfnmet->eigfn = (EIGENFN **)calloc(eigfnmet->nf, sizeof(EIGENFN *));
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        size_t cnum = eigfnmet->eigv[iw].n;
        eigfnmet->eigfn[iw] = (EIGENFN *)calloc(cnum, sizeof(EIGENFN));
    }

    // FFT需要的参数，包括零频和区域外的（从原频散文件确定，因此本模块中 -F 的作用仅用于筛选频段）
    real_t fft_df = eigmet->freqs[1] - eigmet->freqs[0];
    size_t fft_nf = eigmet->nf + 1; // + 1 是为零频
    size_t fft_offset = round(eigmet->freqs[0] / fft_df);
    size_t fft_nt = 2*(fft_nf - 1);
    real_t fft_dt = 1.0 / fft_df;

    // 计算的面波格林函数
    pcplxChnlGrid *grn = (pcplxChnlGrid *) calloc(Ctrl->R.nr, sizeof(*grn));
    pcplxChnlGrid *grn_uiz = (Ctrl->e.active)? (pcplxChnlGrid *) calloc(Ctrl->R.nr, sizeof(*grn_uiz)) : NULL;
    pcplxChnlGrid *grn_uir = (Ctrl->e.active)? (pcplxChnlGrid *) calloc(Ctrl->R.nr, sizeof(*grn_uir)) : NULL;
    for(size_t ir=0; ir<Ctrl->R.nr; ++ir){
        GRT_LOOP_ChnlGrid(im, c){
            grn[ir][im][c] = (cplx_t*)calloc(fft_nf, sizeof(cplx_t));
            if(grn_uiz)  grn_uiz[ir][im][c] = (cplx_t*)calloc(fft_nf, sizeof(cplx_t));
            if(grn_uir)  grn_uir[ir][im][c] = (cplx_t*)calloc(fft_nf, sizeof(cplx_t));
        }
    }

    // Rayleigh or Love
    const int ncols = (eigmet->wtype == GRT_DISPERSION_RAYL)? GRT_RAYL_DIM : GRT_LOVE_DIM;

    // 循环所需要的频率
    #pragma omp parallel for schedule(guided) default(shared)
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        GRT_MODEL1D *local_mod1d = grt_copy_mod1d(mod1d);

        real_t omega = PI2*eigfnmet->freqs[iw];
        EIGENV *eigv = eigfnmet->eigv + iw;

        local_mod1d->omega = omega;
        
        // 假设这些本征值都是从 iref 层的久期函数算出来的，
        // 1. 先根据 iref（z_i+）的垂直波函数，计算出每个物理分界面上侧（z_j-）和下侧（z_j+）的垂直波函数
        // 2. 循环每个深度采样点，使用上侧（z_{j+1}-）和下侧（z_j+）的垂直波函数推导采样点上的垂直波函数以及本征函数
        for(size_t ic = 0; ic < eigv->n; ++ic){
            real_t cphase = eigv->c_roots[ic];
            
            local_mod1d->c_phase = cphase;
            local_mod1d->k = creal(omega)/cphase;
            grt_mod1d_xa_xb(local_mod1d, local_mod1d->k);

            size_t iref = eigv->c_roots_iref[ic];

            // 模型每个物理层介质下方 z_j+ 的垂直波函数
            cplx_t (*mod_potRaylLove_Down)[ncols] = (cplx_t (*)[ncols])calloc(local_mod1d->n, sizeof(cplx_t)*ncols);
            // 模型每个物理层介质上方 z_j- 的垂直波函数，申请 n+1 的内存，方便后续不必讨论“半空间中的上行波场”
            cplx_t (*mod_potRaylLove_Up)[ncols] = (cplx_t (*)[ncols])calloc(local_mod1d->n+1, sizeof(cplx_t)*ncols);

            cplx_t potRaylLove[ncols];  memset(potRaylLove, 0, sizeof(cplx_t)*ncols);
            cplx_t secRaylLove = 0.0;
            grt_secular_function_potential(local_mod1d, cphase, iref, eigfnmet->wtype, &secRaylLove, potRaylLove);
 
            // 计算出每个物理分界面上侧（z_{j+1}-）和 下侧（z_j+）的垂直波函数
            grt_get_mod_potential_Up_Down(local_mod1d, eigfnmet->wtype, ncols, iref, potRaylLove, mod_potRaylLove_Down, mod_potRaylLove_Up);

            // 计算能量积分
            grt_energy_integrals(
                local_mod1d, eigfnmet->wtype, ncols, 
                mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet, &eigfnmet->eigfn[iw][ic]);

            // 计算面波格林函数
            grt_modsum_grn(
                local_mod1d, eigfnmet->wtype, ncols, 
                mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet->eigfn[iw][ic].egyint,
                Ctrl->R.nr, Ctrl->R.rs, iw + fft_offset, 
                Ctrl->e.active, grn, grn_uiz, grn_uir);
            
            GRT_SAFE_FREE_PTR(mod_potRaylLove_Down);
            GRT_SAFE_FREE_PTR(mod_potRaylLove_Up);
        }
        grt_free_mod1d(local_mod1d);
    }


    // 使用fftw3做反傅里叶变换，并保存到 SAC 
    // 其中考虑了升采样倍数
    GRT_FFTW_HOLDER *fh = grt_create_fftw_holder_C2R_1D(fft_nt*Ctrl->W.upsample_n, fft_dt/Ctrl->W.upsample_n, fft_nf, fft_df);

    real_t (* travtPS)[2] = (real_t (*)[2])calloc(Ctrl->R.nr, sizeof(real_t)*2);

    char *s_output_dirprefx = NULL;
    GRT_SAFE_ASPRINTF(&s_output_dirprefx, "%s/%s_%s_%s", Ctrl->O.s_output_dir, modelpath, Ctrl->D.s_depsrc, Ctrl->D.s_deprcv);
    
    // 建立SAC头文件，包含必要的头变量
    SACTRACE *sac = grt_new_SACTRACE(fh->dt, fh->nt, 0.0);
    // 发震时刻作为参考时刻
    sac->hd.o = 0.0; 
    sac->hd.iztype = IO; 
    // 记录震源和台站深度
    sac->hd.evdp = Ctrl->D.depsrc; // km
    sac->hd.stel = (-1.0)*Ctrl->D.deprcv*1e3; // m
    // 写入虚频率
    sac->hd.user0 = 0.0;
    // 写入接受点的Vp,Vs,rho
    sac->hd.user1 = mod1d->Va[mod1d->ircv];
    sac->hd.user2 = mod1d->Vb[mod1d->ircv];
    sac->hd.user3 = mod1d->Rho[mod1d->ircv];
    sac->hd.user4 = mod1d->Qainv[mod1d->ircv];
    sac->hd.user5 = mod1d->Qbinv[mod1d->ircv];
    // 写入震源点的Vp,Vs,rho
    sac->hd.user6 = mod1d->Va[mod1d->isrc];
    sac->hd.user7 = mod1d->Vb[mod1d->isrc];
    sac->hd.user8 = mod1d->Rho[mod1d->isrc];
    
    // 做反傅里叶变换，保存SAC文件
    for(size_t ir = 0; ir < Ctrl->R.nr; ++ir){
        real_t dist = Ctrl->R.rs[ir];
        sac->hd.dist = dist;

        // 文件保存子目录
        char *s_output_subdir = NULL;
        
        GRT_SAFE_ASPRINTF(&s_output_subdir, "%s_%s", s_output_dirprefx, Ctrl->R.s_rs[ir]);
        GRTCheckMakeDir(s_output_subdir);

        // 计算理论走时
        sac->hd.t0 = grt_compute_travt1d(mod1d->Thk, mod1d->Va, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
        strcpy(sac->hd.kt0, "P");
        sac->hd.t1 = grt_compute_travt1d(mod1d->Thk, mod1d->Vb, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
        strcpy(sac->hd.kt1, "S");

        // 时间延迟
        {
            real_t delayT = 0.0;
            if (! Ctrl->E.refFirstP){
                delayT = Ctrl->E.delayT0;
                if(Ctrl->E.delayV0 > 0.0)   delayT += sqrt( GRT_SQUARE(dist) + GRT_SQUARE(Ctrl->D.deprcv - Ctrl->D.depsrc) ) / Ctrl->E.delayV0;
            } else {
                delayT = Ctrl->E.delayT0 + sac->hd.t0;
            }
            // 修改SAC头段时间变量
            sac->hd.b = delayT;
        }

        GRT_LOOP_ChnlGrid(im, c){
            if(! Ctrl->G.doEX  && im==GRT_SRC_M_EX_INDEX)  continue;
            if(! Ctrl->G.doVF  && im==GRT_SRC_M_VF_INDEX)  continue;
            if(! Ctrl->G.doHF  && im==GRT_SRC_M_HF_INDEX)  continue;
            if(! Ctrl->G.doDC  && im>=GRT_SRC_M_DD_INDEX)  continue;

            int modr = GRT_SRC_M_ORDERS[im];
            int sgn=1;  // 用于反转Z分量

            if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;  // 跳过输出0阶的T分量

            // Z分量反转
            sgn = (GRT_ZRT_CODES[c]=='Z') ? -1 : 1;

            char ch = GRT_ZRT_CODES[c];
            
            if(eigfnmet->wtype == GRT_DISPERSION_RAYL && ch =='T')  continue;
            if(eigfnmet->wtype == GRT_DISPERSION_LOVE && ch !='T')  continue;

            write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, 0.0, sac, s_output_subdir, "", sgn, false, grn[ir][im][c]);

            if(Ctrl->e.active){
                write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, 0.0, sac, s_output_subdir, "z", sgn*(-1), false, grn_uiz[ir][im][c]);
                write_one_to_sac(GRT_SRC_M_NAME_ABBR[im], ch, fh, 0.0, sac, s_output_subdir, "r", sgn     , false, grn_uir[ir][im][c]);
            }
        }

        GRT_SAFE_FREE_PTR(s_output_subdir);
        
        // 记录走时
        if(travtPS != NULL){
            travtPS[ir][0] = sac->hd.t0;
            travtPS[ir][1] = sac->hd.t1;
        }
    }

    // 输出警告：当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    if(mod1d->isLiquid[mod1d->isrc]){
        GRTRaiseWarning(
            "The source is located in the liquid layer, "
            "therefore only the Green's Funtions for the Explosion source will be computed.\n");
    }

    grt_free_SACTRACE(sac);

    
    GRT_SAFE_FREE_PTR(modelpath);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
