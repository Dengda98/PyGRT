/**
 * @file   grt_eigenfn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-07
 * 
 *     根据计算的相速度频散(本征值)，求解本征函数（包括位移和应力）
 * 
 */
#include "grt.h"

#include "grt/common/matrix.h"
#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/util.h"
#include "grt/common/search.h"

#include "grt/modal/modal_util.h"
#include "grt/modal/secular.h"
#include "grt/modal/eigenfn.h"
#include "grt/modal/energy.h"


/** 该子模块的参数控制结构体 */
typedef struct {
    /* 输入频散文件 */
    struct {
        bool active;
        char *s_filepath;
        bool readin_period;
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
        real_t *freqs;
        size_t nf;
        bool def_range;  ///< 定义了搜索范围而不是离散的一些点
        bool io_period;  ///< 以周期的形式输入输出
    } F;
    /* 群速度 */
    struct {
        bool active;
        char *s_gdisp_filepath;
    } U;
    /* 敏感核 */
    struct {
        bool active;
        char *s_cpar_filepath;   ///< 相速度对模型参数的偏导数（敏感核）
        char *s_upar_filepath;   ///< 群速度对模型参数的偏导数（敏感核）
        char *s_egyint_filepath;  ///< 输出能量积分
        real_t cpar_dz;   ///< 计算敏感核时对模型深度进行采样
    } K;
    /* 深度采样本征函数 */
    struct {
        bool active;
        real_t *zs;
        size_t nz;
    } Z;
    /* 本征函数输出路径 */
    struct {
        bool active;
        char *s_egn_filepath;
    } O;

    bool calc_egyint;   ///< 是否需要计算能量积分

} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    // C
    GRT_SAFE_FREE_PTR(Ctrl->C.s_filepath);
    // N
    GRT_SAFE_FREE_PTR(Ctrl->N.modes);
    // F
    GRT_SAFE_FREE_PTR(Ctrl->F.freqs);
    // U
    GRT_SAFE_FREE_PTR(Ctrl->U.s_gdisp_filepath);
    // K
    GRT_SAFE_FREE_PTR(Ctrl->K.s_cpar_filepath);
    GRT_SAFE_FREE_PTR(Ctrl->K.s_upar_filepath);
    GRT_SAFE_FREE_PTR(Ctrl->K.s_egyint_filepath);
    // Z
    GRT_SAFE_FREE_PTR(Ctrl->Z.zs);
    // O
    GRT_SAFE_FREE_PTR(Ctrl->O.s_egn_filepath);

    // GRT_SAFE_FREE_PTR_ARRAY(Ctrl->cdisp, Ctrl->F.nf);
    // GRT_SAFE_FREE_PTR_ARRAY(Ctrl->cdisp_iref, Ctrl->F.nf);
    // GRT_SAFE_FREE_PTR(Ctrl->cdisp_num);

    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
    printf("[WAITING TO FINISH]\n");
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, ":C:N::F:U:K:Z:O:sh" )) != -1) {
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

            // 频率值 -Ff1[/f2][/df][+p]  可使用
            case 'F':
                Ctrl->F.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    real_t a1, a2, df;
                    a1 = a2 = 0;
                    df = 1;
                    int nscan = sscanf(token, "%lf/%lf/%lf", &a1, &a2, &df);
                    if( !(nscan == 1 || nscan == 2 || nscan == 3)){
                        GRTBadOptionError(F, "");
                    };
                    if(nscan == 1 && a1 <= 0.0){
                        GRTBadOptionError(F, "Can't set a single nonpositive value.");
                    }
                    if(nscan == 1){
                        a2 = a1;
                        df = 1;
                    }
                    else if(nscan >= 2){
                        if(df <= 0){
                            GRTBadOptionError(F, "Can't set nonpositive df(%lf).", df);
                        }
                        if(a1 < 0.0 || a2 <= 0.0){
                            GRTBadOptionError(F, "Can't set nonpositive f1(%lf), f2(%lf).", a1, a2);
                        }
                        if(a1 > a2){
                            GRTBadOptionError(F, "f1(%lf) > f2(%lf).", a1, a2);
                        }

                        // 跳过零频
                        if(nscan == 3){
                            a1 = GRT_MAX(a1, df);
                            a2 = GRT_MAX(a1, a2);
                        }
                    }

                    bool isperiod = false;
                    // 处理 + 号指令
                    token = strtok(NULL, "+");
                    if(token != NULL){
                        switch (token[0]){
                            case 'p':
                                isperiod = true;
                                break;
                            default:
                                GRTBadOptionError(F, "+%s is not supported.", token);
                                break;
                        }
                    }

                    if(nscan == 3 || nscan == 1){
                        Ctrl->F.nf = floor((a2-a1)/df) + 1;
                        Ctrl->F.freqs = (real_t*)calloc(Ctrl->F.nf, sizeof(real_t));
                        for(size_t i=0; i<Ctrl->F.nf; ++i){
                            if(isperiod){
                                Ctrl->F.freqs[Ctrl->F.nf-1 - i] = 1.0/(a1 + df*i);
                            } else {
                                Ctrl->F.freqs[i] = a1 + df*i;
                            }
                        }
                    }
                    else if(nscan == 2){
                        // 只约束最大值和最小值，内存申请由读取频散文件的过程中进行
                        Ctrl->F.nf = 2;
                        Ctrl->F.freqs = (real_t*)calloc(Ctrl->F.nf, sizeof(real_t));
                        Ctrl->F.freqs[0] = a1;
                        Ctrl->F.freqs[1] = a2;
                        Ctrl->F.def_range = true;
                    }
                    
                    Ctrl->F.io_period = isperiod;

                    GRT_SAFE_FREE_PTR(string);
                }
                break;
            
            // 输出群速度 -U<path>
            case 'U':
                Ctrl->U.active = true;
                Ctrl->U.s_gdisp_filepath = strdup(optarg);
                break;

            // 输出能量积分和敏感核 -K+c<path>+u<path>+z<dz>+x<path>
            case 'K':
                Ctrl->K.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    while(token != NULL){
                        switch (token[0]){
                            case 'c':
                                Ctrl->K.s_cpar_filepath = strdup(token+1);
                                if(strlen(token)==1){
                                    GRTBadOptionError(K, "+%s wrong usage.", token);
                                }
                                break;
                            case 'g':
                                Ctrl->K.s_upar_filepath = strdup(token+1);
                                if(strlen(token)==1){
                                    GRTBadOptionError(K, "+%s wrong usage.", token);
                                }
                                break;                     
                            case 'z':
                                if(1 != sscanf(token+1, "%lf", &Ctrl->K.cpar_dz)){
                                    GRTBadOptionError(K, "+%s wrong usage.", token);
                                }
                                if(Ctrl->K.cpar_dz <= 0.0){
                                    GRTBadOptionError(K, "+%s can't set nonpositive dz(%lf).", token, Ctrl->K.cpar_dz);
                                }
                                break;
                            case 'x':
                                Ctrl->K.s_egyint_filepath = strdup(token+1);
                                if(strlen(token)==1){
                                    GRTBadOptionError(K, "+%s wrong usage.", token);
                                }
                                break;
                            default:
                                GRTBadOptionError(K, "+%s is not supported.", token);
                                break;
                        }
                        token = strtok(NULL, "+");
                    }
                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            // 深度采样  -Zz1[/z2/dz]
            case 'Z':
                Ctrl->Z.active = true;
                {   
                    real_t a1, a2, dif;
                    a1 = a2 = dif = 0;
                    int nscan = sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &dif);
                    if( !(nscan == 1 || nscan == 3)){
                        GRTBadOptionError(Z, "");
                    };
                    if(nscan == 1 && a1 <= 0.0){
                        GRTBadOptionError(Z, "Can't set a single nonpositive value.");
                    }
                    if(nscan == 1){
                        a2 = a1;
                        dif = 1;
                    }
                    else if(nscan == 3){
                        if(dif <= 0){
                            GRTBadOptionError(Z, "Can't set nonpositive dz(%lf).", dif);
                        }
                        if(a1 < 0.0 || a2 <= 0.0){
                            GRTBadOptionError(Z, "Can't set nonpositive z1(%lf), z2(%lf)", a1, a2);
                        }
                        if(a1 > a2){
                            GRTBadOptionError(Z, "z1(%lf) > z2(%lf).", a1, a2);
                        }
                    }

                    Ctrl->Z.nz = floor((a2-a1)/dif) + 1;
                    Ctrl->Z.zs = (real_t*)calloc(Ctrl->Z.nz, sizeof(real_t));
                    for(size_t i=0; i<Ctrl->Z.nz; ++i){
                        Ctrl->Z.zs[i] = a1 + dif*i;
                    }
                }
                break;

            // 输出深度采样的本征函数
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_egn_filepath = strdup(optarg);
                break;

            GRT_Common_Options_in_Switch((char)(optopt));

        }
    }


    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, C);

    // -Z 必须和 -O 同时使用
    if(Ctrl->Z.active + Ctrl->O.active == 1){
        GRTRaiseError("-Z and -O must be set simultaneously.");
    }

    // 是否需要计算能量积分
    Ctrl->calc_egyint = Ctrl->U.active || Ctrl->K.active;

    // -N 的默认项，仅处理基阶
    if( ! Ctrl->N.active ){
        Ctrl->N.nmode = 1;
        Ctrl->N.modes = (size_t *)calloc(Ctrl->N.nmode, sizeof(size_t));
        Ctrl->N.modes[0] = 0;
    }

}



/** 
 * 根据能量积分计算 Rayleight 波和 Love 波的群速度 
 * 
 * @param[in]       omega       圆频率
 * @param[in]       eigenK      本征波数
 * @param[in]       EgyInt      能量积分
 * @param[in]       isRayl      Rayleigh or Love
 * 
 * @return     群速度
 * 
 */
static real_t grt_get_group_velocity(const real_t omega, const real_t eigenK, const cplx_t EgyInt[GRT_EGYINTS_MAX], const DISPER_TYPE wtype)
{
    cplx_t root_c = omega/eigenK;
    cplx_t res;
    if(wtype == GRT_DISPERSION_RAYL){
        res = (EgyInt[1] - EgyInt[2]/(2*eigenK)) / (root_c * EgyInt[0]);
    }
    else if(wtype == GRT_DISPERSION_LOVE){
        res = EgyInt[1] / (root_c * EgyInt[0]);
    }
    else {
        GRTRaiseError("Wrong execution.");
    }
    return creal(res);
}


/** 根据 grt_eigenfn 的一些参数，计算本征函数、能量积分、群速度、相速度敏感核 */
static void grt_eigenfn_egy_gdisp_phaseK_util(
    GRT_MODEL1D *mod1d, EIGENFN_INFO *eigfnmet, const int ncols, 
    const bool calc_egfn, const bool calc_egyint)
{
    // 循环每个频率
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

            // 计算不同深度的本征函数
            if(calc_egfn){
                grt_get_eigenfn_depths(
                    local_mod1d, eigfnmet->wtype, ncols, 
                    mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet, &eigfnmet->eigfn[iw][ic]);
            }

            // 计算能量积分
            if(calc_egyint){
                grt_energy_integrals(
                    local_mod1d, eigfnmet->wtype, ncols, 
                    mod_potRaylLove_Down, mod_potRaylLove_Up, eigfnmet, &eigfnmet->eigfn[iw][ic]);
                // 计算群速度
                eigfnmet->eigv[iw].g_roots[ic] = grt_get_group_velocity(omega, local_mod1d->k, eigfnmet->eigfn[iw][ic].egyint, eigfnmet->wtype);
            }

            GRT_SAFE_FREE_PTR(mod_potRaylLove_Down);
            GRT_SAFE_FREE_PTR(mod_potRaylLove_Up);
        }
        
        grt_free_mod1d(local_mod1d);
    }

    
}


/* 子模块主函数 */
int eigenfn_main(int argc, char **argv){
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
    GRT_SAFE_FREE_PTR(modelpath);

    // 根据命令行参数确定出所需的部分频散信息
    EIGENFN_INFO *eigfnmet = (EIGENFN_INFO *)calloc(1, sizeof(EIGENFN_INFO));
    grt_filter_eigenfn_info(
        Ctrl->F.nf, Ctrl->F.freqs, Ctrl->F.def_range,
        Ctrl->N.nmode, Ctrl->N.modes, eigmet, eigfnmet);
    
    // 初始化 eigfn
    eigfnmet->eigfn = (EIGENFN **)calloc(eigfnmet->nf, sizeof(EIGENFN *));
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        size_t cnum = eigfnmet->eigv[iw].n;
        eigfnmet->eigfn[iw] = (EIGENFN *)calloc(cnum, sizeof(EIGENFN));
    }

    // 判断每个深度值对应的层位
    eigfnmet->nz = Ctrl->Z.nz;
    eigfnmet->zs = (real_t *)calloc(eigfnmet->nz, sizeof(real_t));
    memcpy(eigfnmet->zs, Ctrl->Z.zs, sizeof(real_t)*eigfnmet->nz);
    eigfnmet->z_irefs = (size_t *)calloc(eigfnmet->nz, sizeof(size_t));
    for(size_t iz = 0; iz < Ctrl->Z.nz; ++iz){
        size_t ziref = 0;
        for(ziref=0; ziref < mod1d->n-1; ++ziref){
            if(mod1d->Dep[ziref+1] > Ctrl->Z.zs[iz])  break; 
        }
        eigfnmet->z_irefs[iz] = ziref;
    }

    // Rayleigh or Love
    const int ncols = (eigmet->wtype == GRT_DISPERSION_RAYL)? GRT_RAYL_DIM : GRT_LOVE_DIM;

    if(Ctrl->calc_egyint){
        // 申请 cpar_nz, cpar_zs, cpar_z_irefs
        eigfnmet->cpar_nz = 0;
        real_t thk, h, dz;
        for(size_t iy = 0; iy < mod1d->n; ++iy){
            thk = mod1d->Thk[iy];
            h=0.0;
            dz = (Ctrl->K.cpar_dz==0.0)? thk : Ctrl->K.cpar_dz;
            // 最后一层会略薄于 dz
            while(h < thk){
                // 避免一些因浮点误差而产生的极细层
                if(fabs(thk - h) < dz*1e-5)  break;

                eigfnmet->cpar_z_irefs = (size_t *)realloc(eigfnmet->cpar_z_irefs, sizeof(size_t)*(eigfnmet->cpar_nz+1));
                eigfnmet->cpar_z_irefs[eigfnmet->cpar_nz] = iy;
                eigfnmet->cpar_zs = (real_t *)realloc(eigfnmet->cpar_zs, sizeof(real_t)*(eigfnmet->cpar_nz+1));
                eigfnmet->cpar_zs[eigfnmet->cpar_nz] = mod1d->Dep[iy] + h;

                h += dz;
                eigfnmet->cpar_nz++;

                // 半空间只算一次
                if(iy == mod1d->n-1)  break;
            }
        }
    }

    // 申请敏感核相关内存，赋值 eigfn 其中变量
    for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
        size_t cnum = eigfnmet->eigv[iw].n;
        for(size_t ic = 0; ic < cnum; ++ic){
            EIGENFN *eigfntmp = &eigfnmet->eigfn[iw][ic];

            eigfntmp->eigenC = eigfnmet->eigv[iw].c_roots[ic];
            if(Ctrl->Z.active){
                eigfntmp->fn = (cplx_t (*)[4])calloc(eigfnmet->nz, sizeof(*eigfntmp->fn));
            }
            
            if(Ctrl->calc_egyint){
                eigfntmp->csens = (cplx_t (*)[GRT_SNSTVTY_MAX])calloc(eigfnmet->cpar_nz, sizeof(*eigfntmp->csens));
                eigfntmp->gsens = (cplx_t (*)[GRT_SNSTVTY_MAX])calloc(eigfnmet->cpar_nz, sizeof(*eigfntmp->gsens));
            }
        }
    }

    // // 临时打印 eigfnmet 一些元素测试，看是否能正确读取
    // for(size_t iw = 0; iw < eigfnmet->nf; ++iw){
    //     size_t cnum = eigfnmet->eigv[iw].n;
    //     real_t freq = eigfnmet->freqs[iw];
    //     for(size_t ic = 0; ic < cnum; ++ic){
    //         size_t mode = eigfnmet->modes[ic];

    //         real_t croot = eigfnmet->eigv[iw].c_roots[ic];
    //         size_t ciref = eigfnmet->eigv[iw].c_roots_iref[ic];

    //         printf("# freq=%f(%zu) mode=%zu %e %zu\n", freq, iw, mode, croot, ciref);
    //     }
    // }
    
    // 计算本征函数、能量积分、群速度、相速度敏感核
    grt_eigenfn_egy_gdisp_phaseK_util(mod1d, eigfnmet, ncols, Ctrl->Z.active, Ctrl->calc_egyint);

    // 打印本征函数结果
    if(Ctrl->Z.active){
        grt_output_eigenfns(Ctrl->O.s_egn_filepath, ncols, eigfnmet);
    }

    // 输出群速度
    if(Ctrl->U.active){
        grt_output_gdisp(Ctrl->U.s_gdisp_filepath, eigfnmet);
    }

    // 输出能量积分
    if(Ctrl->K.s_egyint_filepath!=NULL){
        grt_output_energy_integrals(Ctrl->K.s_egyint_filepath, eigfnmet);
    }

    // 输出相速度敏感核
    if(Ctrl->K.s_cpar_filepath!=NULL){
        grt_output_sensitivity(Ctrl->K.s_cpar_filepath, "c", eigfnmet);
    }
    

    // 计算并输出群速度敏感核
    if(Ctrl->K.s_upar_filepath!=NULL){
        grt_group_sensitivity(eigfnmet);
        grt_output_sensitivity(Ctrl->K.s_upar_filepath, "g", eigfnmet);
    }


    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}


