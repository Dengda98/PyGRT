/**
 * @file   grt_strain.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-28
 * 
 *    根据预先合成的位移空间导数，组合成应力张量（由于有衰减，须在频域内进行）
 * 
 */

#include "grt/common/myfftw.h"

#include "grt/common/attenuation.h"
#include "grt/common/sacio2.h"
#include "grt/common/const.h"

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    char *s_synpath;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->s_synpath);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt stress] %s\n\n", GRT_VERSION);printf(
"    Conbine spatial derivatives of displacements into stress tensor.\n"
"    (unit: dyne/cm^2 = 0.1 Pa)\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt stress <syn_dir>\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    (void)Ctrl;
    int opt;
    while ((opt = getopt(argc, argv, ":h")) != -1) {
        switch (opt) {
            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(argc > 1);
}


int stress_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);
    
    // 合成地震图目录路径
    Ctrl->s_synpath = strdup(argv[1]);

    // 检查是否存在该目录
    GRTCheckDirExist(Ctrl->s_synpath);

    // ----------------------------------------------------------------------------------
    // 开始读取计算，输出6个量
    char c1, c2;
    char *s_filepath = NULL;

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = false;
    // 三分量
    const char *chs = NULL;

    // 判断标志性文件是否存在，来判断输出使用ZNE还是ZRT
    GRT_SAFE_ASPRINTF(&s_filepath, "%s/nN.sac", Ctrl->s_synpath);
    rot2ZNE = (access(s_filepath, F_OK) == 0);

    // 指示特定的通道名
    chs = (rot2ZNE)? GRT_ZNE_CODES : GRT_ZRT_CODES;


    // 读取一个头段变量，获得基本参数，分配数组内存
    GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->s_synpath, tolower(chs[0]), chs[0]);
    SACTRACE *insac = grt_read_SACTRACE(s_filepath, true);
    int npts = insac->hd.npts;
    float dt = insac->hd.delta;
    float dist = insac->hd.dist;
    float df = 1.0/(npts*dt);
    int nf = npts/2 + 1;
    float va = insac->hd.user1;
    float vb = insac->hd.user2;
    float rho = insac->hd.user3;
    float Qainv = insac->hd.user4;
    float Qbinv = insac->hd.user5;
    if(va <= 0.0 || vb < 0.0 || rho <= 0.0){
        GRTRaiseError("Bad rcv_va, rcv_vb or rcv_rho in \"%s\" header.\n", s_filepath);
    }
    SACTRACE *outsac = grt_copy_SACTRACE(insac, true);
    grt_free_SACTRACE(insac);

    // 申请内存
    // lamda * 体积应变
    fftwf_complex *lam_ukk = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
    // 不同频率的lambda和mu
    fftwf_complex *lams = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
    fftwf_complex *mus = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*nf);
    // 分配FFTW
    GRT_FFTWF_HOLDER *fwd_fftw_holder = grt_create_fftwf_holder_R2C_1D(npts, dt, nf, df);
    GRT_FFTWF_HOLDER *inv_fftw_holder = grt_create_fftwf_holder_C2R_1D(npts, dt, nf, df);
    // 初始化
    memset(lam_ukk, 0, sizeof(fftwf_complex)*nf);
    memset(lams, 0, sizeof(fftwf_complex)*nf);
    memset(mus, 0, sizeof(fftwf_complex)*nf);
    // 计算不同频率下的拉梅系数
    for(int i=0; i<nf; ++i){
        float freq, w;
        freq = (i==0) ? 0.01f : df*i; // 计算衰减因子不能为0频
        w = PI2 * freq;
        fftwf_complex atta, attb;
        atta = grt_attenuation_law(Qainv, w);
        attb = grt_attenuation_law(Qbinv, w);
        // 乘上1e10，转为dyne/(cm^2)
        mus[i] = vb*vb*attb*attb*rho*1e10;
        lams[i] = va*va*atta*atta*rho*1e10 - 2.0*mus[i];
    }

    // ----------------------------------------------------------------------------------
    // 先计算体积应变u_kk = u_11 + u22 + u33 和 lamda的乘积
    for(int i1=0; i1<3; ++i1){
        c1 = chs[i1];

        // 读取数据 u_{k,k}
        GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->s_synpath, tolower(c1), c1);
        insac = grt_read_SACTRACE(s_filepath, false);
        memcpy(fwd_fftw_holder->w_t, insac->data, sizeof(float)*npts);

        // 累加
        fftwf_execute(fwd_fftw_holder->plan);
        for(int i=0; i<nf; ++i)  lam_ukk[i] += fwd_fftw_holder->W_f[i];
    }
    // 加上协变导数
    if(!rot2ZNE){
        GRT_SAFE_ASPRINTF(&s_filepath, "%s/R.sac", Ctrl->s_synpath);
        insac = grt_read_SACTRACE(s_filepath, false);
        memcpy(fwd_fftw_holder->w_t, insac->data, sizeof(float)*npts);
        fftwf_execute(fwd_fftw_holder->plan);
        for(int i=0; i<nf; ++i)  lam_ukk[i] += fwd_fftw_holder->W_f[i]/dist*1e-5;
    }

    // 乘上lambda系数
    for(int i=0; i<nf; ++i)  lam_ukk[i] *= lams[i];

    // 重新初始化
    grt_reset_fftwf_holder_zero(fwd_fftw_holder);
    grt_reset_fftwf_holder_zero(inv_fftw_holder);

    // ----------------------------------------------------------------------------------
    // 循环6个分量
    for(int i1=0; i1<3; ++i1){
        c1 = chs[i1];
        for(int i2=i1; i2<3; ++i2){
            c2 = chs[i2];

            // 读取数据 u_{i,j}
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->s_synpath, tolower(c2), c1);
            insac = grt_read_SACTRACE(s_filepath, false);
            memcpy(fwd_fftw_holder->w_t, insac->data, sizeof(float)*npts);

            // 累加
            fftwf_execute(fwd_fftw_holder->plan);
            for(int i=0; i<nf; ++i)  inv_fftw_holder->W_f[i] += fwd_fftw_holder->W_f[i];

            // 读取数据 u_{j,i}
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->s_synpath, tolower(c1), c2);
            insac = grt_read_SACTRACE(s_filepath, false);
            memcpy(fwd_fftw_holder->w_t, insac->data, sizeof(float)*npts);
            
            // 累加
            fftwf_execute(fwd_fftw_holder->plan);
            for(int i=0; i<nf; ++i)  inv_fftw_holder->W_f[i] = (inv_fftw_holder->W_f[i] + fwd_fftw_holder->W_f[i]) * mus[i];

            // 对于对角线分量，需加上lambda * u_kk
            if(c1 == c2){
                for(int i=0; i<nf; ++i)  inv_fftw_holder->W_f[i] += lam_ukk[i];
            }

            // 特殊情况需加上协变导数，1e-5是因为km->cm
            if(c1=='R' && c2=='T'){
                // 读取数据 u_T
                GRT_SAFE_ASPRINTF(&s_filepath, "%s/T.sac", Ctrl->s_synpath);
                insac = grt_read_SACTRACE(s_filepath, false);
                memcpy(fwd_fftw_holder->w_t, insac->data, sizeof(float)*npts);
                fftwf_execute(fwd_fftw_holder->plan);
                for(int i=0; i<nf; ++i)  inv_fftw_holder->W_f[i] -= mus[i] * fwd_fftw_holder->W_f[i] / dist * 1e-5;
            }
            else if(c1=='T' && c2=='T'){
                // 读取数据 u_R
                GRT_SAFE_ASPRINTF(&s_filepath, "%s/R.sac", Ctrl->s_synpath);
                insac = grt_read_SACTRACE(s_filepath, false);
                memcpy(fwd_fftw_holder->w_t, insac->data, sizeof(float)*npts);
                fftwf_execute(fwd_fftw_holder->plan);
                for(int i=0; i<nf; ++i)  inv_fftw_holder->W_f[i] += 2.0f * mus[i] * fwd_fftw_holder->W_f[i] / dist * 1e-5;
            }
            
            // 保存到SAC
            fftwf_execute(inv_fftw_holder->plan);
            for(int i=0; i<npts; ++i)  inv_fftw_holder->w_t[i] /= npts;
            memcpy(outsac->data, inv_fftw_holder->w_t, sizeof(float)*npts);
            sprintf(outsac->hd.kcmpnm, "%c%c", c1, c2);
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/stress_%c%c.sac", Ctrl->s_synpath, c1, c2);
            grt_write_SACTRACE(s_filepath, outsac);

            // 置零
            grt_reset_fftwf_holder_zero(inv_fftw_holder);
        }
    }


    grt_destroy_fftwf_holder(fwd_fftw_holder);
    grt_destroy_fftwf_holder(inv_fftw_holder);

    GRT_SAFE_FFTW_FREE_PTR(lam_ukk, f);
    GRT_SAFE_FFTW_FREE_PTR(lams, f);
    GRT_SAFE_FFTW_FREE_PTR(mus, f);

    grt_free_SACTRACE(insac);
    grt_free_SACTRACE(outsac);
    GRT_SAFE_FREE_PTR(s_filepath);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}