/**
 * @file   grt_rftn.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-09
 *
 *    计算平面入射 P-SV 波的接收函数
 *
 */

#include "grt.h"

#include <float.h>

#undef I    ///< 取消标准复数单位宏，避免与命令行选项 I 冲突
#define IMAG _Complex_I   ///< 复数单位，避免与标准复数宏 I 冲突

#define GRT_RFTN_N_ZETA        0.8
#define GRT_RFTN_N_UPSAMPLE    1
#define GRT_RFTN_A_ALP         1.0


/** 入射波类型 */
typedef enum {
    GRT_RFTN_INCIDENT_P = 0,    ///< P 波
    GRT_RFTN_INCIDENT_SV,       ///< SV 波
} GRT_RFTN_INCIDENT;


/** 输出结果类型 */
typedef enum {
    GRT_RFTN_OUTPUT_RATIO = 0,  ///< 接收函数比值
    GRT_RFTN_OUTPUT_Z,          ///< 垂向位移响应
    GRT_RFTN_OUTPUT_R,          ///< 径向位移响应
} GRT_RFTN_OUTPUT;


/** 接收函数子模块的参数控制结构体 */
typedef struct {
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;       ///< 模型路径
        MODEL1D *mod1d;          ///< 模型结构体指针
    } M;
    /** 水平射线参数及入射层信息 */
    struct {
        bool active;
        real_t rayp;             ///< 水平射线参数，单位为 s/km
        size_t incident_idx;     ///< 入射层相对于模型底部的逆序编号
        real_t inca;             ///< 入射角，单位为度
    } P;
    /** 入射角及入射层选项 */
    struct {
        bool active;
        real_t inca;             ///< 入射角，单位为度
        size_t idx;              ///< 入射层相对于模型底部的逆序编号
    } I;
    /** 时间序列及频域计算选项 */
    struct {
        bool active;
        size_t nt;               ///< 时间采样点数
        real_t dt;               ///< 时间采样间隔，单位为 s
        real_t zeta;             ///< 虚频系数
        size_t upsample_n;       ///< 时域上采样倍数
        bool keepAllFreq;        ///< 是否保留全部频率点
        bool skipImagComps;      ///< 是否跳过虚频补偿
    } N;
    /** 入射波类型选项 */
    struct {
        bool active;
        GRT_RFTN_INCIDENT incident; ///< 入射波类型
    } T;
    /** 高斯低通滤波选项 */
    struct {
        bool active;
        real_t alp;              ///< 高斯滤波参数，单位为 Hz
    } A;
    /** 结果时间延迟选项 */
    struct {
        bool active;
        real_t delay;            ///< 用户设置的附加延迟，程序内部会自动进行一部分时延，单位为 s
    } E;
    /** 是否输出绝对位移分量 */
    struct {
        bool active;
        bool write_components;   ///< 是否输出垂向和径向分量
    } W;
    /** 输出目录选项 */
    struct {
        bool active;
        char *s_output_dir;      ///< 输出目录路径
    } O;
    /** 是否静默输出 */
    struct {
        bool active;
    } s;
} GRT_RFTN_CTRL;


/**
 * 释放接收函数模块参数控制结构体
 *
 * @param[in,out] Ctrl  接收函数模块参数控制结构体指针
 */
static void free_Ctrl(GRT_RFTN_CTRL *Ctrl)
{
    if(Ctrl == NULL) return;
    // 释放模型输入
    GRT_SAFE_FREE_PTR(Ctrl->M.s_modelpath);
    grt_free_mod1d(Ctrl->M.mod1d);
    // 释放输出目录
    GRT_SAFE_FREE_PTR(Ctrl->O.s_output_dir);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印使用说明 */
static void print_help(void)
{
printf("\n"
"[grt rftn] %s\n\n", GRT_VERSION); printf(
"    Compute receiver functions for a unit plane incident P or SV wave.\n\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt rftn -M<model> (-P<rayp>|-I<inca>[/<idx>]) -T<P|S>\n"
"             -N<nt>/<dt>[+w<zeta>][+n<fac>][+a][+f]\n"
"             [-A<alp>] [-E<delay>] [-W] -O<output_dir> [-s] [-h]\n\n"
"\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -M<model>    Horizontally layered model file.\n"
"\n"
"    -P<rayp>     Horizontal ray parameter in s/km.\n"
"\n"
"    -I<inca>[/<idx>]\n"
"                 Upgoing incidence angle in degrees at the top of\n"
"                 the reverse-<idx> layer; idx=0 is the halfspace.\n"
"\n"
"    -T<P|S>      Incident P or SV wave.\n"
"\n"
"    -N<nt>/<dt>[+w<zeta>][+n<fac>][+a][+f] \n"
"                 <nt>:   number of points. (NOT requires 2^n).\n"
"                 <dt>:   time interval (secs). \n"
"                 +w<zeta>: define the coefficient of imaginary \n"
"                           frequency wI=zeta*PI/T, where T=nt*dt.\n"
"                           Default zeta=%.1f.\n", GRT_RFTN_N_ZETA); printf(
"                 +n<fac>:  upsampling factor (integer)\n"
"                           i.e.  nt <-- nt * <fac>\n"
"                                 dt <-- dt / <fac>\n"
"                           and calculated frequencies stay unchanged.\n"
"                 +a:       All frequencies are calculated regardless of\n"
"                           how low the frequency is.\n"
"                 +f:       skip the amplitude compensation from \n"
"                           imaginary frequency.\n"
"\n"
"    -A<alp>      Gaussian low-pass filter parameter (default %.1f).\n", GRT_RFTN_A_ALP); printf(
"                 H(f)=exp[-(pi*f/alp)^2]; filter corner frequency ≈ alp/pi.\n"
"\n"
"    -E<delay>    Additional delay before the P arrival. The program\n"
"                 already applies an internal delay; set this only when\n"
"                 more delay is needed.\n"
"\n"
"    -W           Also write absolute Z and R component time series.\n"
"\n"
"    -O<dir>      Output directory.\n"
"\n"
"    -s           Silence progress information.\n"
"\n"
"    -h           Display this help message.\n\n");
}


/** 解析命令行选项 */
static void getopt_from_command(GRT_RFTN_CTRL *Ctrl, int argc, char **argv)
{
    Ctrl->N.zeta = GRT_RFTN_N_ZETA;
    Ctrl->N.upsample_n = GRT_RFTN_N_UPSAMPLE;
    Ctrl->A.alp = GRT_RFTN_A_ALP;

    int opt;
    while((opt = getopt(argc, argv, ":M:P:I:N:T:A:E:O:Wsh")) != -1){
        switch(opt){
            // -M<model>
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                break;

            // -P<rayp>
            case 'P':
                Ctrl->P.active = true;
                {
                    char extra;
                    if(sscanf(optarg, "%lf%c", &Ctrl->P.rayp, &extra) != 1 ||
                        !isfinite(Ctrl->P.rayp) || Ctrl->P.rayp <= 0.0){
                        GRTBadOptionError(P, "rayp must be positive and finite.");
                    }
                }
                break;

            // -I<inca>[/<idx>]
            case 'I':
                {
                    real_t idx = 0.0;
                    int nscan = sscanf(optarg, "%lf/%lf", &Ctrl->I.inca, &idx);
                    if(nscan != 1 && nscan != 2){
                        GRTBadOptionError(I, "Use -I<inca>[/<idx>].");
                    }
                    if(!isfinite(Ctrl->I.inca) || Ctrl->I.inca < 0.0 || Ctrl->I.inca >= 90.0){
                        GRTBadOptionError(I, "incidence angle must be in [0, 90) degrees.");
                    }
                    if(nscan == 2 && (!isfinite(idx) || idx < 0.0 || idx != floor(idx) || idx > SIZE_MAX)){
                        GRTBadOptionError(I, "idx must be a nonnegative integer.");
                    }
                    Ctrl->I.active = true;
                    Ctrl->I.idx = (nscan == 2) ? (size_t)idx : 0;
                }
                break;

            // -N<nt>/<dt>[+w<zeta>][+n<fac>][+a][+f]
            case 'N':
                Ctrl->N.active = true;
                {
                    char *string = strdup(optarg);
                    char *token = strtok(string, "+");
                    if(token == NULL || 2 != sscanf(token, "%zu/%lf", &Ctrl->N.nt, &Ctrl->N.dt) ||
                        Ctrl->N.nt == 0 || Ctrl->N.dt <= 0.0 || !isfinite(Ctrl->N.dt)){
                        GRTBadOptionError(N, "");
                    }

                    token = strtok(NULL, "+");
                    while(token != NULL){
                        switch(token[0]){
                            // +n<fac>
                            case 'n':
                                if(1 != sscanf(token + 1, "%zu", &Ctrl->N.upsample_n) || Ctrl->N.upsample_n == 0){
                                    GRTBadOptionError(N, "");
                                }
                                break;

                            // +w<zeta>
                            case 'w':
                                if(1 != sscanf(token + 1, "%lf", &Ctrl->N.zeta) || Ctrl->N.zeta <= 0.0){
                                    GRTBadOptionError(N, "");
                                }
                                break;

                            // +a
                            case 'a':
                                if(token[1] != '\0') GRTBadOptionError(N, "");
                                Ctrl->N.keepAllFreq = true;
                                break;

                            // +f
                            case 'f':
                                if(token[1] != '\0') GRTBadOptionError(N, "");
                                Ctrl->N.skipImagComps = true;
                                break;
                            default:
                                GRTBadOptionError(N, "+%s is not supported.", token);
                        }
                        token = strtok(NULL, "+");
                    }
                    GRT_SAFE_FREE_PTR(string);
                }
                break;

            // -T<P|S>
            case 'T':
                if(strcmp(optarg, "P") == 0){
                    Ctrl->T.incident = GRT_RFTN_INCIDENT_P;
                }
                else if(strcmp(optarg, "S") == 0){
                    Ctrl->T.incident = GRT_RFTN_INCIDENT_SV;
                }
                else{
                    GRTBadOptionError(T, "Use -TP or -TS.");
                }
                Ctrl->T.active = true;
                break;

            // -A<alp>
            case 'A':
                Ctrl->A.active = true;
                {
                    char extra;
                    if(sscanf(optarg, "%lf%c", &Ctrl->A.alp, &extra) != 1 ||
                        !isfinite(Ctrl->A.alp) || Ctrl->A.alp <= 0.0){
                        GRTBadOptionError(A, "alp must be positive and finite.");
                    }
                }
                break;

            // -E<delay>
            case 'E':
                Ctrl->E.active = true;
                {
                    char extra;
                    if(sscanf(optarg, "%lf%c", &Ctrl->E.delay, &extra) != 1 ||
                        !isfinite(Ctrl->E.delay)){
                        GRTBadOptionError(E, "delay must be finite.");
                    }
                }
                break;

            // -O<dir>
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.s_output_dir = strdup(optarg);
                break;

            // -W
            case 'W':
                Ctrl->W.active = true;
                Ctrl->W.write_components = true;
                break;

            // -s
            case 's':
                Ctrl->s.active = true;
                break;
            GRT_Common_Options_in_Switch((char)optopt);
        }
    }

    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, M);
    GRTCheckOptionActive(Ctrl, N);
    GRTCheckOptionActive(Ctrl, T);
    GRTCheckOptionActive(Ctrl, O);
    if(Ctrl->P.active == Ctrl->I.active){
        GRTRaiseError("Set exactly one of -P<rayp> and -I<inca>[/<idx>].");
    }
    if(optind != argc){
        GRTRaiseError("Unexpected positional argument %s. Use '-h' for help.", argv[optind]);
    }
}


/**
 * 计算自由表面位移响应矩阵
 *
 * R_EV = (D11 + D12 RU) (I - RD_RL RU)^-1 TU_RL
 *
 * 行对应 q_m 和 w_m，列对应单位振幅的 P 波和 SV 波入射
 *
 * @param[in,out] mstat  模型状态结构体
 * @param[out]    R_EV   自由表面位移响应矩阵
 * @return              矩阵求逆状态
 */
static int compute_surface_REV(MODEL1D_STATE *mstat, cplx_t R_EV[2][2])
{
    RT_MATRIX M_RL = {0};
    RT_MATRIX M_interface = {0};
    RT_MATRIX M_next = {0};
    grt_reset_RT_matrix_PSV(&M_RL);

    // 从底部向上递推所有物理层的整体 R/T 矩阵
    for(size_t iy = 1; iy < mstat->mod1d->n; ++iy){
        grt_reset_RT_matrix_PSV(&M_interface);
        grt_RT_matrix_PSV(mstat, iy, &M_interface);
        if(M_interface.stats == GRT_INVERSE_FAILURE) return GRT_INVERSE_FAILURE;

        grt_delay_RT_matrix_PSV(mstat, iy, &M_interface);

        grt_reset_RT_matrix_PSV(&M_next);
        grt_recursion_RT_matrix_PSV(&M_RL, &M_interface, &M_next);
        if(M_next.stats == GRT_INVERSE_FAILURE) return GRT_INVERSE_FAILURE;

        M_RL = M_next;
    }

    grt_topbound_RU_PSV(mstat);
    if(mstat->M_top.stats == GRT_INVERSE_FAILURE){
        return GRT_INVERSE_FAILURE;
    }

    // 计算自由表面到顶层向上传播波的位移转换矩阵
    cplx_t D11[2][2], D12[2][2], R_EV_surface[2][2];
    grt_get_layer_D11(mstat, 0, D11);
    grt_get_layer_D12(mstat, 0, D12);
    grt_cmat2x2_mul(D12, mstat->M_top.RU, R_EV_surface);
    grt_cmat2x2_add(D11, R_EV_surface, R_EV_surface);

    // 求解底部反射波和顶界面反射波之间的多次反射反馈
    cplx_t feedback[2][2], inv_feedback[2][2], upgoing[2][2];
    grt_cmat2x2_mul(M_RL.RD, mstat->M_top.RU, feedback);
    grt_cmat2x2_one_sub(feedback);
    if(grt_cmat2x2_inv(feedback, inv_feedback) == GRT_INVERSE_FAILURE){
        return GRT_INVERSE_FAILURE;
    }
    grt_cmat2x2_mul(inv_feedback, M_RL.TU, upgoing);
    grt_cmat2x2_mul(R_EV_surface, upgoing, R_EV);

    return GRT_INVERSE_SUCCESS;
}


/**
 * 截取入射波能够传播到的模型部分
 *
 * @param[in] mod1d   原始模型结构体指针
 * @param[in] nlayer  截取后的模型层数
 * @return            截取后的模型结构体指针
 */
static MODEL1D *truncate_model(MODEL1D *mod1d, const size_t nlayer)
{
    if(nlayer == mod1d->n) return mod1d;
    if(nlayer == 0 || nlayer > mod1d->nmodarr){
        GRTRaiseError("Invalid incident-layer truncation (%zu layers).", nlayer);
    }

    MODEL1D *truncated = grt_read_mod1d_from_modarr(
        nlayer, (const real_t (*)[GRT_MODARR_NCOL])mod1d->modarr, -1.0, -1.0, true);
    if(truncated == NULL){
        GRTRaiseError("Failed to truncate the model at the incident layer.");
    }
    grt_free_mod1d(mod1d);
    return truncated;
}


/**
 * 计算指定波型在模型中的垂向传播时间
 *
 * @param[in] mod1d     模型结构体指针
 * @param[in] rayp      水平射线参数
 * @param[in] is_p_wave 是否计算 P 波传播时间
 * @return              垂向传播时间
 */
static real_t compute_vertical_travel_time(const MODEL1D *mod1d, const real_t rayp, const bool is_p_wave)
{
    real_t travel_time = 0.0;
    for(size_t i = 0; i + 1 < mod1d->n; ++i){
        real_t velocity = is_p_wave ? mod1d->Va[i] : mod1d->Vb[i];
        if(velocity <= 0.0){
            GRTRaiseError(
                "Cannot compute a %s travel time through layer %zu with nonpositive velocity.",
                is_p_wave ? "P-wave" : "S-wave", i + 1);
        }
        real_t vertical_slowness_square = 1.0/(velocity*velocity) - rayp*rayp;
        if(vertical_slowness_square < 0.0){
            GRTRaiseError(
                "The %s ray is not propagating in layer %zu for rayp %.9g.",
                is_p_wave ? "P" : "S", i + 1, rayp);
        }
        travel_time += mod1d->Thk[i] * sqrt(vertical_slowness_square);
    }
    return travel_time;
}


/** 根据水平射线参数或入射角准备入射层和截取后的模型 */
static void prepare_incidence(GRT_RFTN_CTRL *Ctrl)
{
    MODEL1D *mod1d = Ctrl->M.mod1d;
    size_t original_nlayer = mod1d->n;
    size_t keep_nlayer;

    if(Ctrl->I.active){
        if(Ctrl->I.idx >= original_nlayer){
            GRTRaiseError("Incident-layer idx %zu is out of range [0, %zu].", Ctrl->I.idx, original_nlayer - 1);
        }
        keep_nlayer = original_nlayer - Ctrl->I.idx;
        Ctrl->P.incident_idx = Ctrl->I.idx;

        size_t incident_layer = keep_nlayer - 1;
        real_t velocity = (Ctrl->T.incident == GRT_RFTN_INCIDENT_P) ?
            mod1d->Va[incident_layer] : mod1d->Vb[incident_layer];
        if(velocity <= 0.0){
            GRTRaiseError("The selected incident layer has no valid %s-wave velocity.",
                Ctrl->T.incident == GRT_RFTN_INCIDENT_P ? "P" : "SV");
        }
        Ctrl->P.rayp = sin(Ctrl->I.inca*DEG1) / velocity;
    }
    else{
        keep_nlayer = 0;
        for(size_t i = 0; i < mod1d->n; ++i){
            if(Ctrl->P.rayp*mod1d->Va[i] < 1.0){
                keep_nlayer = i + 1;
            }
            else{
                break;
            }
        }
        if(keep_nlayer == 0){
            GRTRaiseError("No real P-wave path exists for rayp %.9g.", Ctrl->P.rayp);
        }
        Ctrl->P.incident_idx = original_nlayer - keep_nlayer;
    }

    Ctrl->M.mod1d = truncate_model(mod1d, keep_nlayer);

    real_t incident_velocity = (Ctrl->T.incident == GRT_RFTN_INCIDENT_P) ?
        Ctrl->M.mod1d->Va[Ctrl->M.mod1d->n - 1] : Ctrl->M.mod1d->Vb[Ctrl->M.mod1d->n - 1];
    real_t sin_inca = Ctrl->P.rayp * incident_velocity;
    if(sin_inca >= 1.0){
        GRTRaiseError(
            "The incident %s wave is not propagating in the selected layer for rayp %.9g.",
            Ctrl->T.incident == GRT_RFTN_INCIDENT_P ? "P" : "SV", Ctrl->P.rayp);
    }
    Ctrl->P.inca = asin(sin_inca) / DEG1;

}


/**
 * 计算一个复频率下的自由表面位移响应
 *
 * @param[in,out] mstat             模型状态结构体
 * @param[in]     k0                复水平波数
 * @param[in]     incident          入射波类型
 * @param[out]    q_m               自由表面水平位移势响应
 * @param[out]    w_m               自由表面垂向位移势响应
 * @param[out]    incident_velocity 入射波速度因子
 * @return                          矩阵求逆状态
 */
static int compute_one_frequency(
    MODEL1D_STATE *mstat, cplx_t k0, GRT_RFTN_INCIDENT incident,
    cplx_t *q_m, cplx_t *w_m, cplx_t *incident_velocity)
{
    if(!isfinite(creal(k0)) || !isfinite(cimag(k0)) || cabs(k0) == 0.0){
        return GRT_INVERSE_FAILURE;
    }

    grt_update_mod1d_state_k(mstat, k0);

    cplx_t R_EV[2][2] = {{0}};
    if(compute_surface_REV(mstat, R_EV) == GRT_INVERSE_FAILURE){
        return GRT_INVERSE_FAILURE;
    }

    size_t incident_index = (incident == GRT_RFTN_INCIDENT_P) ? 0 : 1;
    *q_m = R_EV[0][incident_index];
    *w_m = R_EV[1][incident_index];
    if(!isfinite(creal(*q_m)) || !isfinite(cimag(*q_m)) ||
        !isfinite(creal(*w_m)) || !isfinite(cimag(*w_m))){
        return GRT_INVERSE_FAILURE;
    }

    size_t incident_layer = mstat->mod1d->n - 1;
    if(incident == GRT_RFTN_INCIDENT_P){
        *incident_velocity = mstat->mod1d->Va[incident_layer] * mstat->atna[incident_layer];
    }
    else{
        *incident_velocity = mstat->mod1d->Vb[incident_layer] * mstat->atnb[incident_layer];
    }

    return GRT_INVERSE_SUCCESS;
}


/**
 * 在单次频率循环中记录自由表面响应
 *
 * @param[in,out] mod1d             模型结构体指针
 * @param[in]     rayp              水平射线参数
 * @param[in]     incident          入射波类型
 * @param[out]    q_m               自由表面水平位移势响应数组
 * @param[out]    w_m               自由表面垂向位移势响应数组
 * @param[out]    incident_velocity 入射波速度因子数组
 * @param[in]     nf                频率点数
 * @param[in]     df                频率间隔
 * @param[in]     wI                虚频系数
 */
static void compute_frequency_responses(
    MODEL1D *mod1d, const real_t rayp, const GRT_RFTN_INCIDENT incident,
    cplx_t q_m[], cplx_t w_m[], cplx_t incident_velocity[],
    const size_t nf, const real_t df, const real_t wI)
{
    mod1d->omgref = PI2*(nf - 1)*df;
    const cplx_t invalid = NAN + IMAG*NAN;

    for(size_t iw = 0; iw < nf; ++iw){
        q_m[iw] = invalid;
        w_m[iw] = invalid;
        incident_velocity[iw] = invalid;
    }

    for(size_t iw = 0; iw < nf; ++iw){
        real_t freq = iw*df;
        cplx_t omega = PI2*freq - IMAG*wI;
        cplx_t k0 = omega*rayp;

        MODEL1D_STATE *mstat = grt_init_mod1d_state(mod1d);
        grt_update_mod1d_state_omega(mstat, omega, false);
        int status = compute_one_frequency(
            mstat, k0, incident, &q_m[iw], &w_m[iw], &incident_velocity[iw]);
        grt_free_mod1d_state(mstat);
        if(status == GRT_INVERSE_FAILURE) continue;
    }
}


/**
 * 根据自由表面响应组装一种输出类型的频谱
 *
 * @param[in]  q_m               自由表面水平位移势响应数组
 * @param[in]  w_m               自由表面垂向位移势响应数组
 * @param[in]  incident_velocity 入射波速度因子数组
 * @param[in]  incident          入射波类型
 * @param[in]  output            输出结果类型
 * @param[out] spectrum          频率响应数组
 * @param[in]  nf                频率点数
 * @param[in]  df                频率间隔
 * @param[in]  wI                虚频系数
 * @param[in]  alp               高斯滤波参数
 * @return                       组装失败的频率点数
 */
static size_t assemble_spectrum(
    const cplx_t q_m[], const cplx_t w_m[], const cplx_t incident_velocity[],
    const GRT_RFTN_INCIDENT incident, const GRT_RFTN_OUTPUT output,
    cplx_t spectrum[], const size_t nf, const real_t df,
    const real_t wI, const real_t alp)
{
    memset(spectrum, 0, nf*sizeof(*spectrum));
    size_t nfailed = 0;
    const real_t gauss_exponent_limit = -log(DBL_MIN);

    for(size_t iw = 0; iw < nf; ++iw){
        if(!isfinite(creal(q_m[iw])) || !isfinite(cimag(q_m[iw])) ||
            !isfinite(creal(w_m[iw])) || !isfinite(cimag(w_m[iw]))){
            ++nfailed;
            continue;
        }

        real_t gauss_exponent = GRT_SQUARE(PI*(iw*df)/alp);
        if(gauss_exponent > gauss_exponent_limit) continue;
        real_t gauss = exp(-gauss_exponent);

        cplx_t response;
        if(output == GRT_RFTN_OUTPUT_RATIO){
            cplx_t denominator = (incident == GRT_RFTN_INCIDENT_P) ? w_m[iw] : q_m[iw];
            if(cabs(denominator) == 0.0){
                ++nfailed;
                continue;
            }
            response = (incident == GRT_RFTN_INCIDENT_P) ? IMAG*q_m[iw]/w_m[iw] : -IMAG*w_m[iw]/q_m[iw];
        }
        else{
            cplx_t omega = PI2*(iw*df) - IMAG*wI;
            if(!isfinite(creal(incident_velocity[iw])) ||
                !isfinite(cimag(incident_velocity[iw])) || cabs(omega) == 0.0){
                ++nfailed;
                continue;
            }

            // 将单位入射波势系数转换为单位位移
            cplx_t unit_displacement_factor = incident_velocity[iw]/(IMAG*omega);
            cplx_t displacement = (output == GRT_RFTN_OUTPUT_Z) ? w_m[iw] : IMAG*q_m[iw];
            if(incident == GRT_RFTN_INCIDENT_SV){
                // 应用 SV 入射波相位修正
                displacement *= IMAG;
            }
            response = unit_displacement_factor * displacement;
        }

        spectrum[iw] = gauss*response;
        if(!isfinite(creal(spectrum[iw])) || !isfinite(cimag(spectrum[iw]))){
            spectrum[iw] = 0.0;
            ++nfailed;
        }
    }

    return nfailed;
}


/**
 * 逆变换频谱并写入 SAC 文件
 *
 * @param[in] incident          入射波类型
 * @param[in] output            输出结果类型
 * @param[in] path              SAC 文件路径
 * @param[in] spectrum          频率响应数组
 * @param[in] nt                时间采样点数
 * @param[in] dt                时间采样间隔
 * @param[in] rayp              水平射线参数
 * @param[in] alp               高斯滤波参数
 * @param[in] skip_imag_comps   是否跳过虚频补偿
 * @param[in] nf                频率点数
 * @param[in] df                频率间隔
 * @param[in] wI                虚频系数
 * @param[in] begin_time        SAC 时间序列起始时间
 * @param[in] phase_time        频谱相位延迟
 */
static void write_result_sac(
    const GRT_RFTN_INCIDENT incident, const GRT_RFTN_OUTPUT output, const char *path,
    const cplx_t spectrum[], const size_t nt, const real_t dt,
    const real_t rayp, const real_t alp, const bool skip_imag_comps,
    const size_t nf, const real_t df, const real_t wI,
    const real_t begin_time, const real_t phase_time)
{
    GRT_FFTW_HOLDER *fh = grt_create_fftw_holder_C2R_1D(nt, dt, nf, df);
    SACTRACE *sac = grt_new_SACTRACE(fh->dt, (int)fh->nt, begin_time);

    sac->hd.o = 0.0;
    sac->hd.iztype = IO;
    sac->hd.user0 = wI;
    sac->hd.user1 = rayp;
    sac->hd.user2 = alp;
    sac->hd.user3 = (float)incident;
    sac->hd.user4 = (float)output;
    snprintf(sac->hd.kuser0, sizeof(sac->hd.kuser0), "wI");
    snprintf(sac->hd.kuser1, sizeof(sac->hd.kuser1), "rayp");
    snprintf(sac->hd.kuser2, sizeof(sac->hd.kuser2), "alp");
    const char incident_name = (incident == GRT_RFTN_INCIDENT_P) ? 'P' : 'S';
    if(output == GRT_RFTN_OUTPUT_RATIO){
        snprintf(sac->hd.kcmpnm, sizeof(sac->hd.kcmpnm), "%c_RFTN", incident_name);
    }
    else if(output == GRT_RFTN_OUTPUT_Z){
        snprintf(sac->hd.kcmpnm, sizeof(sac->hd.kcmpnm), "%c_Z", incident_name);
    }
    else{
        snprintf(sac->hd.kcmpnm, sizeof(sac->hd.kcmpnm), "%c_R", incident_name);
    }

    grt_reset_fftw_holder_zero(fh);
    // 将相位延迟乘到频谱上，使时间序列从指定时刻开始
    cplx_t phase = exp(IMAG*PI2*df*phase_time);
    cplx_t phase_factor = 1.0;
    for(size_t iw = 0; iw < nf; ++iw){
        fh->W_f[iw] = spectrum[iw]*phase_factor;
        phase_factor *= phase;
    }

    if(fh->naive_inv){
        grt_naive_inverse_transform_double(fh);
    }
    else {
        fftw_execute(fh->plan);
    }

    // 仅在未跳过虚频补偿时恢复复频率变换的指数因子
    real_t coefficient = df;
    real_t imag_factor = 1.0;
    if(!skip_imag_comps){
        coefficient *= exp(phase_time*wI);
        imag_factor = exp(wI*fh->dt);
    }
    for(size_t it = 0; it < fh->nt; ++it){
        sac->data[it] = (float)(fh->w_t[it]*coefficient);
        coefficient *= imag_factor;
    }

    if(grt_write_SACTRACE(path, sac) != 0){
        GRTRaiseError("Failed to write SAC file %s.", path);
    }

    grt_free_SACTRACE(sac);
    grt_destroy_fftw_holder(fh);
}


/** rftn 子模块主函数 */
int rftn_main(int argc, char **argv)
{
    GRT_RFTN_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    Ctrl->M.mod1d = grt_read_mod1d_from_file(Ctrl->M.s_modelpath, -1.0, -1.0, true);
    if(Ctrl->M.mod1d == NULL){
        free_Ctrl(Ctrl);
        return EXIT_FAILURE;
    }
    prepare_incidence(Ctrl);
    if(Ctrl->T.incident == GRT_RFTN_INCIDENT_SV && Ctrl->M.mod1d->isLiquid[Ctrl->M.mod1d->n - 1]){
        GRTRaiseError("SV incidence is not defined when the bottom halfspace is liquid.");
    }

    real_t winT = Ctrl->N.nt*Ctrl->N.dt;
    if(!isfinite(winT) || winT <= 0.0){
        GRTRaiseError("The time window must be finite and positive.");
    }
    size_t nf = Ctrl->N.nt/2 + 1;
    real_t df = 1.0/winT;
    real_t wI = Ctrl->N.zeta*PI/winT;
    // 高斯滤波函数的逆变换与 exp[-(alp*t)^2] 成正比
    real_t delay = Ctrl->E.delay - 2.0/Ctrl->A.alp;
    real_t tp;
    real_t ts = NAN;
    real_t begin_time;
    real_t phase_components;
    if(Ctrl->T.incident == GRT_RFTN_INCIDENT_P){
        tp = compute_vertical_travel_time(Ctrl->M.mod1d, Ctrl->P.rayp, true);
        begin_time = delay;
    }
    else{
        ts = compute_vertical_travel_time(Ctrl->M.mod1d, Ctrl->P.rayp, false);
        tp = compute_vertical_travel_time(Ctrl->M.mod1d, Ctrl->P.rayp, true);
        begin_time = tp - ts + delay;
    }
    phase_components = tp + delay;

    if(!Ctrl->s.active){
        const char *wave_name = (Ctrl->T.incident == GRT_RFTN_INCIDENT_P) ? "P" : "SV";
        const char *layer_name = (Ctrl->P.incident_idx == 0) ? "halfspace" : "layer";

        GRTRaiseInfo("Incident wave: %s", wave_name);
        GRTRaiseInfo("rayp = %.9g s/km", Ctrl->P.rayp);
        GRTRaiseInfo("incidence layer: reverse-%zu %s", Ctrl->P.incident_idx, layer_name);
        GRTRaiseInfo("incidence angle = %.6g deg.", Ctrl->P.inca);

        if(Ctrl->T.incident == GRT_RFTN_INCIDENT_P){
            GRTRaiseInfo("Vertical travel time: Tp = %.6g s", tp);
            GRTRaiseInfo("begin time = %.6g s", begin_time);
        }
        else{
            GRTRaiseInfo("Vertical travel times: Tp = %.6g s", tp);
            GRTRaiseInfo("Vertical travel times: Ts = %.6g s", ts);
            GRTRaiseInfo("begin time = %.6g s", begin_time);
        }
    }
    GRTCheckMakeDir(Ctrl->O.s_output_dir);

    size_t nt = Ctrl->N.nt*Ctrl->N.upsample_n;
    real_t dt = Ctrl->N.dt/Ctrl->N.upsample_n;
    cplx_t *q_m = calloc(nf, sizeof(*q_m));
    cplx_t *w_m = calloc(nf, sizeof(*w_m));
    cplx_t *incident_velocity = calloc(nf, sizeof(*incident_velocity));
    cplx_t *spectrum = calloc(nf, sizeof(*spectrum));

    // 频率递推只执行一次，后续输出复用 q_m 和 w_m
    compute_frequency_responses(
        Ctrl->M.mod1d, Ctrl->P.rayp, Ctrl->T.incident,
        q_m, w_m, incident_velocity, nf, df, wI);

    const GRT_RFTN_OUTPUT output_types[] = {
        GRT_RFTN_OUTPUT_RATIO, GRT_RFTN_OUTPUT_Z, GRT_RFTN_OUTPUT_R,
    };
    const char *output_names[] = {"rftn", "Z", "R"};
    const real_t phase_times[] = {begin_time, phase_components, phase_components};
    size_t noutputs = Ctrl->W.write_components ? 3 : 1;
    const char incident_name = (Ctrl->T.incident == GRT_RFTN_INCIDENT_P) ? 'P' : 'S';
    for(size_t ioutput = 0; ioutput < noutputs; ++ioutput){
        GRT_RFTN_OUTPUT output = output_types[ioutput];
        size_t nfailed = assemble_spectrum(
            q_m, w_m, incident_velocity, Ctrl->T.incident, output,
            spectrum, nf, df, wI, Ctrl->A.alp);

        if(nfailed > 0 && !Ctrl->s.active){
            GRTRaiseWarning("%zu frequency points were singular or non-finite and were set to zero.", nfailed);
        }

        char *path = NULL;
        GRT_SAFE_ASPRINTF(&path, "%s/%c_%s.sac", Ctrl->O.s_output_dir, incident_name, output_names[ioutput]);
        write_result_sac(
            Ctrl->T.incident, output, path, spectrum, nt, dt,
            Ctrl->P.rayp, Ctrl->A.alp, Ctrl->N.skipImagComps,
            nf, df, wI, begin_time, phase_times[ioutput]);

        if(!Ctrl->s.active){
            if(output == GRT_RFTN_OUTPUT_RATIO){
                GRTRaiseInfo("Saved receiver function to %s", path);
            }
            else if(output == GRT_RFTN_OUTPUT_Z){
                GRTRaiseInfo("Saved vertical response to %s", path);
            }
            else{
                GRTRaiseInfo("Saved radial response to %s", path);
            }
        }
        GRT_SAFE_FREE_PTR(path);
    }

    GRT_SAFE_FREE_PTR(q_m);
    GRT_SAFE_FREE_PTR(w_m);
    GRT_SAFE_FREE_PTR(incident_velocity);
    GRT_SAFE_FREE_PTR(spectrum);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
