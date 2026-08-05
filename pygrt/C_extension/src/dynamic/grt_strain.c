/**
 * @file   grt_strain.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-28
 * 
 *    根据预先合成的位移空间导数，组合成应变张量
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
"[grt strain] %s\n\n", GRT_VERSION);printf(
"    Conbine spatial derivatives of displacements into strain tensor.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt strain <syn_dir>\n"
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

void grt_compute_strain(
    size_t npts, float dist, float *const u[GRT_CHANNEL_NUM],
    float *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    float *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE)
{
    const char *chs = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;

    for(size_t i=0; i<npts; ++i){
        // 联络项（1e-5: km→cm）：r≠0 用 u/r；r=0 改用 ∂_r u，与 syn 中 (1/r)∂_θ 有限部分配套
        float ur_over_r = GRT_IS_ZERO(dist) ? upar[1][1][i] : (u[1][i] / dist * 1e-5f);
        float ut_over_r = GRT_IS_ZERO(dist) ? upar[1][2][i] : (u[2][i] / dist * 1e-5f);

        for(int c=0; c<GRT_CHANNEL_NUM; ++c){
            for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
                float val = 0.5f * (upar[c2][c][i] + upar[c][c2][i]);
                if(chs[c]=='R' && chs[c2]=='T'){
                    val -= 0.5f * ut_over_r;
                }
                else if(chs[c]=='T' && chs[c2]=='T'){
                    val += ur_over_r;
                }
                res[c2][c][i] = val;
            }
        }
    }
}



/** 子模块主函数 */
int strain_main(int argc, char **argv){
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
    float dist = insac->hd.dist;
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
    grt_compute_strain(npts, dist, u, upar, res, rot2ZNE);

    // 写出6个分量
    for(int i1=0; i1<3; ++i1){
        c1 = chs[i1];
        for(int i2=i1; i2<3; ++i2){
            c2 = chs[i2];
            memcpy(outsac->data, res[i2][i1], sizeof(*outsac->data)*npts);
            sprintf(outsac->hd.kcmpnm, "%c%c", c1, c2);
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/strain_%c%c.sac", Ctrl->s_synpath, c1, c2);
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
