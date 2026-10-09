/**
 * @file   grt_coulomb.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-03
 *
 * 动态库仑应力：tau_s + friction * sigma_n
 */
#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    struct {
        bool active;
        char *path;  ///< 输入结果目录
    } G;
    struct {
        bool active;
        real_t friction;  ///< 有效摩擦系数
    } F;
} GRT_MODULE_CTRL;

/** 释放参数控制结构体 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->G.path);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(void)
{
printf("\n"
"[grt coulomb] %s\n\n", GRT_VERSION);printf(
"    Compute dynamic Coulomb stress change from sigma_n.sac and tau_s.sac,\n"
"    which are normally produced by `sproj`. The result is written to each\n"
"    receiver directory as coulomb.sac = tau_s + friction * sigma_n\n"
"    (unit: dyne/cm^2 = 0.1 Pa). Existing results are overwritten.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt coulomb -G<directory> -F<friction> [-h]\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<directory> Input receiver directory or root of multiple receiver\n"
"                  directories containing sigma_n.sac and tau_s.sac.\n"
"\n"
"    -F<friction>  Nonnegative dimensionless effective friction coefficient.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
);
}

/**
 * 读取并检查命令行参数
 *
 * @param[out] Ctrl  参数控制结构体
 * @param[in]  argc  命令行参数数量
 * @param[in]  argv  命令行参数数组
 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    // 解析输入目录和有效摩擦系数
    int opt;
    while((opt = getopt(argc, argv, ":G:F:h")) != -1) {
        switch(opt) {
            case 'G':
                Ctrl->G.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->G.path);
                Ctrl->G.path = strdup(optarg);
                break;

            case 'F': {
                char extra;
                if(sscanf(optarg, "%lf%c", &Ctrl->F.friction, &extra) != 1 || Ctrl->F.friction < 0) {
                    GRTRaiseError("Friction must be nonnegative.");
                }
                Ctrl->F.active = true;
                break;
            }
            GRT_Common_Options_in_Switch((char)optopt);
        }
    }

    // 输入投影结果目录和有效摩擦系数均为必选项
    GRTCheckOptionActive(Ctrl, G);
    if(optind != argc) {
        GRTRaiseError("Unexpected positional argument %s. Use '-h' for help.", argv[optind]);
    }
    GRTCheckOptionActive(Ctrl, F);
}

/** 子模块主函数 */
int coulomb_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    // 将单台目录或多台根目录展开为接收目录列表
    size_t nr;
    char **dirs = grt_syn_output_receiver_directories(Ctrl->G.path, &nr);

    // 各接收目录独立读取投影结果并计算库仑应力
    for(size_t ir = 0; ir < nr; ++ir) {
        char *path = NULL;
        GRT_SAFE_ASPRINTF(&path, "%s/sigma_n.sac", dirs[ir]);
        SACTRACE *normal = grt_read_SACTRACE(path, false);
        GRT_SAFE_FREE_PTR(path);
        GRT_SAFE_ASPRINTF(&path, "%s/tau_s.sac", dirs[ir]);
        SACTRACE *shear = grt_read_SACTRACE(path, false);
        GRT_SAFE_FREE_PTR(path);

        // 只有使用同一时间网格的两条波形才能逐样本组合
        if(normal->hd.npts <= 0 || normal->hd.delta <= 0 || 
           normal->hd.npts != shear->hd.npts || normal->hd.delta != shear->hd.delta || normal->hd.b != shear->hd.b) {
            GRTRaiseError("Invalid or inconsistent stress time sampling at %s.", dirs[ir]);
        }

        // 逐时刻将摩擦加权的法向应力加到滑动方向的剪应力上
        for(int n = 0; n < shear->hd.npts; ++n) {
            shear->data[n] += Ctrl->F.friction * normal->data[n];
        }

        // 继承剪应力的接收头段，更新分量名并保存到当前接收目录
        snprintf(shear->hd.kcmpnm, sizeof(shear->hd.kcmpnm), "coulomb");
        GRT_SAFE_ASPRINTF(&path, "%s/coulomb.sac", dirs[ir]);
        grt_write_SACTRACE(path, shear);
        GRT_SAFE_FREE_PTR(path);
        grt_free_SACTRACE(normal);
        grt_free_SACTRACE(shear);
        GRT_SAFE_FREE_PTR(dirs[ir]);
    }
    GRT_SAFE_FREE_PTR(dirs);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
