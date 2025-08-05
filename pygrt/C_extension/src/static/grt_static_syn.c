/**
 * @file   grt_static_syn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-02-18
 * 
 *    根据计算好的静态格林函数，定义震源机制以及方位角等，生成合成的静态三分量位移场
 * 
 */


#include "common/const.h"
#include "common/radiation.h"
#include "common/coord.h"

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    char *name;
    /** 旋转到 Z, N, E */
    struct {
        bool active;
    } N;
    /** 放大系数 */
    struct {
        bool active;
        bool mult_src_mu;
        MYREAL M0;
        MYREAL src_mu;
    } S;  
    /** 剪切源 */
    struct {
        bool active;
    } M;
    /** 单力源 */
    struct {
        bool active;
    } F;
    /** 矩张量源 */
    struct {
        bool active;
    } T;
    /** 是否计算空间导数 */
    struct {
        bool active;
    } e;

    // 存储不同震源的震源机制相关参数的数组
    MYREAL mchn[MECHANISM_NUM];

    // 方向因子数组
    MYREAL srcRadi[SRC_M_NUM][CHANNEL_NUM];

    // 最终要计算的震源类型
    MYINT computeType;
    char s_computeType[3];

} GRT_MODULE_CTRL;

/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    free(Ctrl->name);
    free(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt static syn] %s\n\n", GRT_VERSION);printf(
"    Compute static displacement with the outputs of \n"
"    module `static_greenfn` (reading from stdin).\n"
"    Three components are:\n"
"       + Up (Z),\n"
"       + Radial Outward (R),\n"
"       + Transverse Clockwise (T),\n"
"    and the units are cm. You can add -N to rotate ZRT to ZNE.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static syn -S[u]<scale> \n"
"              [-M<strike>/<dip>/<rake>]\n"
"              [-T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>]\n"
"              [-F<fn>/<fe>/<fz>] \n"
"              [-N] [-e]\n"
"              < <grn>\n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -S[u]<scale>  Scale factor to all kinds of source. \n"
"                  + For Explosion, Shear and Moment Tensor,\n"
"                    unit of <scale> is dyne-cm. \n"
"                  + For Single Force, unit of <scale> is dyne.\n"
"                  + Since \"\\mu\" exists in scalar seismic moment\n"
"                    (\\mu*A*D), you can simply set -Su<scale>, <scale>\n"
"                    equals A*D (Area*Slip, [cm^3]), and <scale> will \n"
"                    multiply \\mu automatically in program.\n"
"\n"
"    For source type, you can only set at most one of\n"
"    '-M', '-T' and '-F'. If none, an Explosion is used.\n"
"\n"
"    -M<strike>/<dip>/<rake>\n"
"                  Three angles to define a fault. \n"
"                  The angles are in degree.\n"
"\n"
"    -T<Mxx>/<Mxy>/<Mxz>/<Myy>/<Myz>/<Mzz>\n"
"                  Six elements of Moment Tensor. \n"
"                  x (North), y (East), z (Downward).\n"
"                  Notice they will be scaled by <scale>.\n"
"\n"
"    -F<fn>/<fe>/<fz>\n"
"                  North, East and Vertical(Downward) Forces.\n"
"                  Notice they will be scaled by <scale>.\n"
"\n"
"    -N            Components of results will be Z, N, E.\n"
"\n"
"    -e            Compute the spatial derivatives, ui_z and ui_r,\n"
"                  of displacement u. In filenames, prefix \"r\" means \n"
"                  ui_r and \"z\" means ui_z. \n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    Say you have computed Static Green's functions with following command:\n"
"        grt static greenfn -Mmilrow -D2/0 -X-5/5/10 -Y-5/5/10 > grn\n"
"\n"
"    Then you can get static displacement of Explosion\n"
"        grt static syn -Su1e16 < grn > syn_ex\n"
"\n"
"    or Shear\n"
"        grt static syn -Su1e16 -M100/20/80 < grn > syn_dc\n"
"\n"
"    or Single Force\n"
"        grt static syn -S1e20 -F0.5/-1.2/3.3 < grn > syn_sf\n"
"\n"
"    or Moment Tensor\n"
"        grt static syn -Su1e16 -T2.3/0.2/-4.0/0.3/0.5/1.2 < grn > syn_mt\n"
"\n\n\n"
"\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    const char *command = Ctrl->name;

    // 先为个别参数设置非0初始值
    Ctrl->computeType = GRT_SYN_COMPUTE_EX;
    sprintf(Ctrl->s_computeType, "%s", "EX");

    int opt;
    while ((opt = getopt(argc, argv, ":S:M:F:T:Neh")) != -1) {
        switch (opt) {
            // 放大系数
            case 'S':
                Ctrl->S.active = true;
                {   
                    // 检查是否存在字符u，若存在表明需要乘上震源处的剪切模量
                    char *upos=NULL;
                    if((upos=strchr(optarg, 'u')) != NULL){
                        Ctrl->S.mult_src_mu = true;
                        *upos = ' ';
                    }
                }
                if(0 == sscanf(optarg, "%lf", &Ctrl->S.M0)){
                    GRTBadOptionError(command, S, "");
                };
                break;

            // 剪切震源
            case 'M':
                Ctrl->M.active = true;
                Ctrl->computeType = GRT_SYN_COMPUTE_DC;
                {
                    double strike, dip, rake;
                    sprintf(Ctrl->s_computeType, "%s", "DC");
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &strike, &dip, &rake)){
                        GRTBadOptionError(command, M, "");
                    };
                    if(strike < 0.0 || strike > 360.0){
                        GRTBadOptionError(command, M, "Strike must be in [0, 360].");
                    }
                    if(dip < 0.0 || dip > 90.0){
                        GRTBadOptionError(command, M, "Dip must be in [0, 90].");
                    }
                    if(rake < -180.0 || rake > 180.0){
                        GRTBadOptionError(command, M, "Rake must be in [-180, 180].");
                    }
                    Ctrl->mchn[0] = strike;
                    Ctrl->mchn[1] = dip;
                    Ctrl->mchn[2] = rake;
                }
                break;

            // 单力源
            case 'F':
                Ctrl->F.active = true;
                Ctrl->computeType = GRT_SYN_COMPUTE_SF;
                {
                    double fn, fe, fz;
                    sprintf(Ctrl->s_computeType, "%s", "SF");
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &fn, &fe, &fz)){
                        GRTBadOptionError(command, F, "");
                    };
                    Ctrl->mchn[0] = fn;
                    Ctrl->mchn[1] = fe;
                    Ctrl->mchn[2] = fz;
                }
                break;

            // 张量震源
            case 'T':
                Ctrl->T.active = true;
                Ctrl->computeType = GRT_SYN_COMPUTE_MT;
                {
                    double Mxx, Mxy, Mxz, Myy, Myz, Mzz;
                    sprintf(Ctrl->s_computeType, "%s", "MT");
                    if(6 != sscanf(optarg, "%lf/%lf/%lf/%lf/%lf/%lf", &Mxx, &Mxy, &Mxz, &Myy, &Myz, &Mzz)){
                        GRTBadOptionError(command, T, "");
                    };
                    Ctrl->mchn[0] = Mxx;
                    Ctrl->mchn[1] = Mxy;
                    Ctrl->mchn[2] = Mxz;
                    Ctrl->mchn[3] = Myy;
                    Ctrl->mchn[4] = Myz;
                    Ctrl->mchn[5] = Mzz;
                }
                break;

            // 是否计算位移空间导数, 影响 calcUTypes 变量
            case 'e':
                Ctrl->e.active = true;
                break;

            // 是否旋转到ZNE, 影响 rot2ZNE 变量
            case 'N':
                Ctrl->N.active = true;
                break;

            GRT_Common_Options_in_Switch(command, (char)(optopt));
        }
    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(command, argc > 1);
    GRTCheckOptionActive(command, Ctrl, S);

    // 只能使用一种震源
    if(Ctrl->M.active + Ctrl->F.active + Ctrl->T.active > 1){
        GRTRaiseError("[%s] Error! Only support at most one of \"-M\", \"-F\" and \"-T\". Use \"-h\" for help.\n", command);
    }
}




/** 子模块主函数 */
int static_syn_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);
    const char *command = Ctrl->name;

    getopt_from_command(Ctrl, argc, argv);

    // 辐射因子
    // double srcRadi[SRC_M_NUM][CHANNEL_NUM]={0};

    // 从标准输入中读取静态格林函数表
    double x0, y0, grn[SRC_M_NUM][CHANNEL_NUM]={0}, syn[CHANNEL_NUM]={0}, syn_upar[CHANNEL_NUM][CHANNEL_NUM]={0};
    double grn_uiz[SRC_M_NUM][CHANNEL_NUM]={0}, grn_uir[SRC_M_NUM][CHANNEL_NUM]={0};

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = Ctrl->N.active;

    // 根据参数设置，选择分量名
    const char *chs = (rot2ZNE)? ZNEchs : ZRTchs;


    // 建立一个指针数组，方便读取多列数据
    const int max_ncol = 47;
    double *pt_grn[max_ncol];
    // 按照特定顺序
    {
        double **pt = &pt_grn[0];
        *(pt++) = &x0;
        *(pt++) = &y0;
        for(int m=0; m<3; ++m){
            for(int k=0; k<SRC_M_NUM; ++k){
                for(int c=0; c<CHANNEL_NUM; ++c){
                    if(SRC_M_ORDERS[k]==0 && ZRTchs[c] == 'T')  continue;

                    if(m==0){
                        *pt = &grn[k][c];
                    } else if(m==1){
                        *pt = &grn_uiz[k][c];
                    } else if(m==2){
                        *pt = &grn_uir[k][c];
                    }
                    pt++;
                }
            }
        }
    }

    // 是否已打印输出的列名
    bool printHead = false;

    // 输入列数
    int ncols = 0;

    // 方位角
    double azrad = 0.0;

    // 物性参数
    double src_va=0.0, src_vb=0.0, src_rho=0.0, src_mu=0.0;
    double rcv_va=0.0, rcv_vb=0.0, rcv_rho=0.0;

    // 用于计算位移空间导数的比例系数
    double upar_scale=1.0; 

    // 计算和位移相关量的种类（1-位移，2-ui_z，3-ui_r，4-ui_t）
    int calcUTypes = (Ctrl->e.active)? 4 : 1;

    // 震中距
    double dist = 0.0;

    char line[1024];
    int iline = 0;
    while(fgets(line, sizeof(line), stdin)){
        iline++;
        if(iline == 1){
            // 读取震源物性参数
            if(3 != sscanf(line, "# %lf %lf %lf", &src_va, &src_vb, &src_rho)){
                GRTRaiseError("[%s] Error! Unable to read src property from \"%s\". \n", command, line);
            }
            if(src_va <= 0.0 || src_vb < 0.0 || src_rho <= 0.0){
                GRTRaiseError("[%s] Error! Bad src_va, src_vb or src_rho from \"%s\". \n", command, line);
            }
            if(src_vb == 0.0 && Ctrl->S.mult_src_mu){
                GRTRaiseError("[%s] Error! Zero src_vb from \"%s\". "
                    "Maybe you try to use -Su<scale> but the source is in the liquid. "
                    "Use -S<scale> instead.\n", command, line);
                exit(EXIT_FAILURE);
            }
            src_mu = src_vb*src_vb*src_rho*1e10;

            if(Ctrl->S.mult_src_mu)  Ctrl->S.M0 *= src_mu;
        }
        else if(iline == 2){
            // 读取场点物性参数
            if(3 != sscanf(line, "# %lf %lf %lf", &rcv_va, &rcv_vb, &rcv_rho)){
                GRTRaiseError("[%s] Error! Unable to read rcv property from \"%s\". \n", command, line);
            }
            if(rcv_va <= 0.0 || rcv_vb < 0.0 || rcv_rho <= 0.0){
                GRTRaiseError("[%s] Error! Bad rcv_va, rcv_vb or rcv_rho in line %d from \"%s\". \n", command, iline, line);
            }
        }
        else if(iline == 3){
            // 根据列长度判断是否有位移空间导数
            char *copyline = strdup(line+1);  // +1去除首个#字符
            char *token = strtok(copyline, " ");
            while (token != NULL) {
                ncols++;
                token = strtok(NULL, " ");
            }
            free(copyline);

            // 想合成位移空间导数但输入的格林函数没有
            if(Ctrl->e.active && ncols < max_ncol){
                GRTRaiseError("[%s] Error! The input has no spatial derivatives. \n", command);
            }
        }
        if(line[0] == '#')  continue;

        // 读取该行数据
        char *copyline = strdup(line);
        char *token = strtok(copyline, " ");
        for(int i=0; i<ncols; ++i){
            sscanf(token, "%lf", pt_grn[i]);  token = strtok(NULL, " ");
        }
        free(copyline);

        // 计算方位角
        azrad = atan2(y0, x0);

        // 计算震中距
        dist = sqrt(x0*x0 + y0*y0);
        if(dist < 1e-5)  dist=1e-5;

        // 先输出列名
        if(!printHead){
            // 打印物性参数
            fprintf(stdout, "# "GRT_REAL_FMT" "GRT_REAL_FMT" "GRT_REAL_FMT"\n", src_va, src_vb, src_rho);
            fprintf(stdout, "# "GRT_REAL_FMT" "GRT_REAL_FMT" "GRT_REAL_FMT"\n", rcv_va, rcv_vb, rcv_rho);
            
            char XX[20];
            sprintf(XX, GRT_STRING_FMT, "X(km)"); XX[0]='#';
            fprintf(stdout, "%s", XX);
            fprintf(stdout, GRT_STRING_FMT, "Y(km)");
            char s_channel[5];
            for(int i=0; i<CHANNEL_NUM; ++i){
                sprintf(s_channel, "%s%c", Ctrl->s_computeType, toupper(chs[i])); 
                fprintf(stdout, GRT_STRING_FMT, s_channel);
            }

            if(Ctrl->e.active){
                for(int k=0; k<CHANNEL_NUM; ++k){
                    for(int i=0; i<CHANNEL_NUM; ++i){
                        sprintf(s_channel, "%c%s%c", tolower(chs[k]), Ctrl->s_computeType, toupper(chs[i])); 
                        fprintf(stdout, GRT_STRING_FMT, s_channel);
                    }
                }
            }

            fprintf(stdout, "\n");
            printHead = true;
        }

        double (*grn3)[CHANNEL_NUM];  // 使用对应类型的格林函数
        double tmpsyn[CHANNEL_NUM];
        for(int ityp=0; ityp<calcUTypes; ++ityp){
            // 求位移空间导数时，需调整比例系数
            switch (ityp){
                // 合成位移
                case 0:
                    upar_scale=1.0;
                    break;

                // 合成ui_z
                case 1:
                // 合成ui_r
                case 2:
                    upar_scale=1e-5;
                    break;

                // 合成ui_t
                case 3:
                    upar_scale=1e-5 / dist;
                    break;
                    
                default:
                    break;
            }

            if(ityp==1){
                grn3 = grn_uiz;
            } else if(ityp==2){
                grn3 = grn_uir;
            } else {
                grn3 = grn;
            }

            tmpsyn[0] = tmpsyn[1] = tmpsyn[2] = 0.0;
            // 计算震源辐射因子
            set_source_radiation(Ctrl->srcRadi, Ctrl->computeType, ityp==3, Ctrl->S.M0, upar_scale, azrad, Ctrl->mchn);

            for(int i=0; i<CHANNEL_NUM; ++i){
                for(int k=0; k<SRC_M_NUM; ++k){
                    tmpsyn[i] += grn3[k][i] * Ctrl->srcRadi[k][i];
                }
            }

            // 保存数据
            for(int i=0; i<CHANNEL_NUM; ++i){
                if(ityp == 0){
                    syn[i] = tmpsyn[i];
                } else {
                    syn_upar[ityp-1][i] = tmpsyn[i];
                }
            }
        }

        // 是否要转到ZNE
        if(rot2ZNE){
            if(Ctrl->e.active){
                rot_zrt2zxy_upar(azrad, syn, syn_upar, dist*1e5);
            } else {
                rot_zxy2zrt_vec(-azrad, syn);
            }
        }

        // 输出数据
        fprintf(stdout, GRT_REAL_FMT GRT_REAL_FMT, x0, y0);
        for(int i=0; i<CHANNEL_NUM; ++i){
            fprintf(stdout, GRT_REAL_FMT, syn[i]);
        }
        if(Ctrl->e.active){
            for(int i=0; i<CHANNEL_NUM; ++i){
                for(int k=0; k<CHANNEL_NUM; ++k){
                    fprintf(stdout, GRT_REAL_FMT, syn_upar[i][k]);
                }
            }
        }
        

        fprintf(stdout, "\n");
        
    }

    if(iline==0){
        GRTRaiseError("[%s] Error! Empty input. \n", command);
    }

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}