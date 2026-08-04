/**
 * @file   grt_strain.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-28
 * 
 *    根据预先合成的位移空间导数，组合成应力张量（由于有衰减，须在频域内进行）
 * 
 */

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

void grt_compute_stress(
    size_t npts, float dt, float dist, float va, float vb, float rho,
    float Qainv, float Qbinv, float *const u[GRT_CHANNEL_NUM],
    float *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    float *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE)
{
    const char *chs = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;
    size_t nf = npts/2 + 1;
    float df = 1.0f/(npts*dt);
    fftwf_complex *lam_ukk = fftwf_malloc(sizeof(*lam_ukk)*nf);
    fftwf_complex *lams = fftwf_malloc(sizeof(*lams)*nf);
    fftwf_complex *mus = fftwf_malloc(sizeof(*mus)*nf);
    GRT_FFTWF_HOLDER *fwd = grt_create_fftwf_holder_R2C_1D(npts, dt, nf, df);
    GRT_FFTWF_HOLDER *inv = grt_create_fftwf_holder_C2R_1D(npts, dt, nf, df);

    memset(lam_ukk, 0, sizeof(*lam_ukk)*nf);
    for(size_t i=0; i<nf; ++i){
        float freq = (i==0) ? 0.01f : df*i;
        float w = PI2 * freq;
        fftwf_complex atta = grt_attenuation_law(Qainv, PI2*(nf-1)*df, w);
        fftwf_complex attb = grt_attenuation_law(Qbinv, PI2*(nf-1)*df, w);
        mus[i] = vb*vb*attb*attb*rho*1e10f;
        lams[i] = va*va*atta*atta*rho*1e10f - 2.0f*mus[i];
    }

    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        memcpy(fwd->w_t, upar[c][c], sizeof(float)*npts);
        fftwf_execute(fwd->plan);
        for(size_t i=0; i<nf; ++i)  lam_ukk[i] += fwd->W_f[i];
    }

    // ZRT 联络项 u/r（1e-5: km→cm）；时域先算好再 FFT，避免频域再分支缩放
    float *ur_over_r = (float *)malloc(sizeof(float)*npts);
    float *ut_over_r = (float *)malloc(sizeof(float)*npts);
    for(size_t i=0; i<npts; ++i){
        ur_over_r[i] = u[1][i] / dist * 1e-5f;
        ut_over_r[i] = u[2][i] / dist * 1e-5f;
    }

    if(!rot2ZNE){
        memcpy(fwd->w_t, ur_over_r, sizeof(float)*npts);
        fftwf_execute(fwd->plan);
        for(size_t i=0; i<nf; ++i)  lam_ukk[i] += fwd->W_f[i];
    }
    for(size_t i=0; i<nf; ++i)  lam_ukk[i] *= lams[i];

    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
            memcpy(fwd->w_t, upar[c2][c], sizeof(float)*npts);
            fftwf_execute(fwd->plan);
            for(size_t i=0; i<nf; ++i)  inv->W_f[i] += fwd->W_f[i];

            memcpy(fwd->w_t, upar[c][c2], sizeof(float)*npts);
            fftwf_execute(fwd->plan);
            for(size_t i=0; i<nf; ++i)  inv->W_f[i] = (inv->W_f[i] + fwd->W_f[i]) * mus[i];
            if(c == c2){
                for(size_t i=0; i<nf; ++i)  inv->W_f[i] += lam_ukk[i];
            }
            if(chs[c]=='R' && chs[c2]=='T'){
                memcpy(fwd->w_t, ut_over_r, sizeof(float)*npts);
                fftwf_execute(fwd->plan);
                for(size_t i=0; i<nf; ++i)  inv->W_f[i] -= mus[i]*fwd->W_f[i];
            }
            else if(chs[c]=='T' && chs[c2]=='T'){
                memcpy(fwd->w_t, ur_over_r, sizeof(float)*npts);
                fftwf_execute(fwd->plan);
                for(size_t i=0; i<nf; ++i)  inv->W_f[i] += 2.0f*mus[i]*fwd->W_f[i];
            }
            fftwf_execute(inv->plan);
            for(size_t i=0; i<npts; ++i)  res[c2][c][i] = inv->w_t[i]/npts;
            grt_reset_fftwf_holder_zero(inv);
        }
    }

    GRT_SAFE_FREE_PTR(ur_over_r);
    GRT_SAFE_FREE_PTR(ut_over_r);

    grt_destroy_fftwf_holder(fwd);
    grt_destroy_fftwf_holder(inv);
    fftwf_free(lam_ukk);
    fftwf_free(lams);
    fftwf_free(mus);
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

    float *u[GRT_CHANNEL_NUM];
    float *upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    float *res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c.sac", Ctrl->s_synpath, chs[c]);
        insac = grt_read_SACTRACE(s_filepath, false);
        u[c] = insac->data;
        insac->data = NULL;
        grt_free_SACTRACE(insac);
        for(int c2=0; c2<GRT_CHANNEL_NUM; ++c2){
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->s_synpath, tolower(chs[c2]), chs[c]);
            insac = grt_read_SACTRACE(s_filepath, false);
            upar[c2][c] = insac->data;
            insac->data = NULL;
            grt_free_SACTRACE(insac);
            res[c2][c] = calloc(npts, sizeof(*res[c2][c]));
        }
    }
    grt_compute_stress(npts, dt, dist, va, vb, rho, Qainv, Qbinv, u, upar, res, rot2ZNE);

    // 写出6个分量
    for(int i1=0; i1<3; ++i1){
        c1 = chs[i1];
        for(int i2=i1; i2<3; ++i2){
            c2 = chs[i2];
            memcpy(outsac->data, res[i2][i1], sizeof(*outsac->data)*npts);
            sprintf(outsac->hd.kcmpnm, "%c%c", c1, c2);
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/stress_%c%c.sac", Ctrl->s_synpath, c1, c2);
            grt_write_SACTRACE(s_filepath, outsac);
        }
    }

    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        free(u[c]);
        for(int c2=0; c2<GRT_CHANNEL_NUM; ++c2){
            free(upar[c2][c]);
            free(res[c2][c]);
        }
    }
    grt_free_SACTRACE(outsac);
    GRT_SAFE_FREE_PTR(s_filepath);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
