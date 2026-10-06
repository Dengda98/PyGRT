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
    struct {
        bool active;
        char *path;  ///< 合成结果目录
    } G;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->G.path);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt strain] %s\n\n", GRT_VERSION);printf(
"    Combine spatial derivatives of displacements into strain tensor.\n"
"    Input must be a `syn` output computed with -e.\n"
"    Six SAC files are written to the same directory.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt strain -G<syn_dir> [-h]\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<syn_dir>   Input receiver directory or root of multiple receiver\n"
"                  directories produced by `syn` with -e. Results are\n"
"                  written to each receiver directory.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
);
}


/** 从命令行中读取选项，处理后记录到参数控制结构体 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    // 仅为兼容旧命令保留模块名后的目录参数，帮助文档统一使用 -G
    bool legacy_path = argc > 1 && argv[1][0] != '-';
    if(legacy_path) {
        Ctrl->G.active = true;
        Ctrl->G.path = strdup(argv[1]);
        optind = 2;
    }

    int opt;
    while((opt = getopt(argc, argv, ":G:h")) != -1) {
        switch(opt) {
            case 'G':
                if(legacy_path) {
                    GRTRaiseError("Cannot combine -G with a positional directory. Use '-h' for help.");
                }
                Ctrl->G.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->G.path);
                Ctrl->G.path = strdup(optarg);
                break;

            GRT_Common_Options_in_Switch((char)optopt);
        }
    }

    GRTCheckOptionActive(Ctrl, G);
    if(optind != argc) {
        GRTRaiseError("Unexpected positional argument %s. Use '-h' for help.", argv[optind]);
    }
    if(legacy_path) {
        GRTRaiseWarning("Passing the input directory as a positional argument is deprecated. Use \"grt %s -G<syn_dir>\" instead.", GRT_MODULE_NAME);
    }
}

/** 由位移偏导合成应变张量 */
static void compute_strain(
    size_t npts, real_t dist, real_t *const u[GRT_CHANNEL_NUM],
    real_t *const upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM],
    real_t *const res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM], bool rot2ZNE)
{
    const char *chs = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;

    for(size_t i=0; i<npts; ++i){
        // 联络项（1e-5: km→cm）：r≠0 用 u/r；r=0 改用 ∂_r u，与 syn 中 (1/r)∂_θ 有限部分配套
        real_t ur_over_r = GRT_IS_ZERO(dist) ? upar[1][1][i] : (u[1][i] / dist * 1e-5);
        real_t ut_over_r = GRT_IS_ZERO(dist) ? upar[1][2][i] : (u[2][i] / dist * 1e-5);

        for(int c=0; c<GRT_CHANNEL_NUM; ++c){
            for(int c2=c; c2<GRT_CHANNEL_NUM; ++c2){
                real_t val = 0.5 * (upar[c2][c][i] + upar[c][c2][i]);
                if(chs[c]=='R' && chs[c2]=='T'){
                    val -= 0.5 * ut_over_r;
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
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    // ----------------------------------------------------------------------------------
    // 开始读取计算，输出6个量
    char c1, c2;
    char *s_filepath = NULL;

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = false;
    // 三分量
    const char *chs = NULL;

    // 判断标志性文件是否存在，来判断输出使用ZNE还是ZRT
    GRT_SAFE_ASPRINTF(&s_filepath, "%s/nN.sac", Ctrl->G.path);
    rot2ZNE = (access(s_filepath, F_OK) == 0);

    // 指示特定的通道名
    chs = (rot2ZNE)? GRT_ZNE_CODES : GRT_ZRT_CODES;

    // 读取一个头段变量，获得基本参数，分配数组内存
    GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->G.path, tolower(chs[0]), chs[0]);
    SACTRACE *insac = grt_read_SACTRACE(s_filepath, true);
    int npts = insac->hd.npts;
    real_t dist = insac->hd.dist;
    SACTRACE *outsac = grt_copy_SACTRACE(insac, true);
    grt_free_SACTRACE(insac);

    real_t *u[GRT_CHANNEL_NUM];
    real_t *upar[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    real_t *res[GRT_CHANNEL_NUM][GRT_CHANNEL_NUM];
    for(int c=0; c<GRT_CHANNEL_NUM; ++c){
        GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c.sac", Ctrl->G.path, chs[c]);
        insac = grt_read_SACTRACE(s_filepath, false);
        u[c] = insac->data;
        insac->data = NULL;
        grt_free_SACTRACE(insac);
        for(int c2=0; c2<GRT_CHANNEL_NUM; ++c2){
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->G.path, tolower(chs[c2]), chs[c]);
            insac = grt_read_SACTRACE(s_filepath, false);
            upar[c2][c] = insac->data;
            insac->data = NULL;
            grt_free_SACTRACE(insac);
            res[c2][c] = GRT_SAFE_CALLOC(npts, sizeof(*res[c2][c]));
        }
    }
    compute_strain(npts, dist, u, upar, res, rot2ZNE);

    // 写出6个分量
    for(int i1=0; i1<3; ++i1){
        c1 = chs[i1];
        for(int i2=i1; i2<3; ++i2){
            c2 = chs[i2];
            memcpy(outsac->data, res[i2][i1], sizeof(*outsac->data)*npts);
            sprintf(outsac->hd.kcmpnm, "%c%c", c1, c2);
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/strain_%c%c.sac", Ctrl->G.path, c1, c2);
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
