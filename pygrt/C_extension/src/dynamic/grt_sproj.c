/**
 * @file   grt_sproj.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-03
 *
 * 动态应力张量在接收断层上的逐时刻投影
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
        char *path;  ///< 接收点文件
    } Q;
    struct {
        bool active;
        bool full;      ///< 是否给出完整接收机制
        bool force;     ///< 是否强制替换已有滑动角
        real_t strike;  ///< 接收走向，度
        real_t dip;     ///< 接收倾角，度
        real_t rake;    ///< 接收滑动角，度
    } M;
} GRT_MODULE_CTRL;

/** 释放参数控制结构体 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->G.path);
    GRT_SAFE_FREE_PTR(Ctrl->Q.path);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(void)
{
printf("\n"
"[grt sproj] %s\n\n", GRT_VERSION);printf(
"    Project dynamic stress tensors onto receiver-fault geometry.\n"
"    The input directory must contain the six stress_*.sac files produced\n"
"    by `stress`. Both ZNE and ZRT stress components are supported.\n"
"    The results are written to each receiver directory as sigma_n.sac\n"
"    and tau_s.sac (unit: dyne/cm^2 = 0.1 Pa). Existing results are\n"
"    overwritten. Normal stress is positive in tension; the positive\n"
"    shear direction is defined by rake.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt sproj -G<directory> [-M<geometry> | -Q<file>] [-h]\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<directory> Input receiver directory or root of multiple receiver\n"
"                  directories.\n"
"\n"
"    -M<geometry>  By default, read receiver geometry from SAC headers.\n"
"                  For ordinary points, use <strike>/<dip>/<rake> in degree.\n"
"                  For finite receivers, use <rake> to fill undefined rakes\n"
"                  or <rake>+f to replace the rake of every receiver.\n"
"\n"
"    -Q<file>      For ordinary points, read replacement receiver geometry.\n"
"                  Each row must contain \"north east depth strike dip rake\".\n"
"                  Coordinates (km) and point count must match SAC headers\n"
"                  in receiver-directory index order, starting from zero.\n"
"                  All points are checked before writing any results.\n"
"                  Mutually exclusive with -M.\n"
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
    int opt;
    while((opt = getopt(argc, argv, ":G:Q:M:h")) != -1) {
        switch(opt) {
            case 'G':
                Ctrl->G.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->G.path);
                Ctrl->G.path = strdup(optarg);
                break;

            case 'Q':
                Ctrl->Q.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->Q.path);
                Ctrl->Q.path = strdup(optarg);
                break;

            case 'M': {
                Ctrl->M.active = true;

                // 将可选的 +f 标记与主参数分开，主参数可以是完整机制或单独的滑动角
                char *value = strdup(optarg), *suffix = strchr(value, '+'), extra;
                if(suffix) {
                    if(strcmp(suffix, "+f")) {
                        GRTRaiseError("Only +f is supported after -M.");
                    }
                    *suffix = '\0';
                    Ctrl->M.force = true;
                }

                int n = sscanf(value, "%lf/%lf/%lf%c", &Ctrl->M.strike, &Ctrl->M.dip, &Ctrl->M.rake, &extra);
                if(n == 1) {
                    if(sscanf(value, "%lf%c", &Ctrl->M.rake, &extra) != 1) {
                        GRTRaiseError("Invalid -M rake.");
                    }
                } else if(n == 3) {
                    Ctrl->M.full = true;
                } else {
                    GRTRaiseError("Expected -Mstrike/dip/rake or -Mrake[+f].");
                }
                GRT_SAFE_FREE_PTR(value);

                // 检查滑动角范围，完整机制还需要检查走向和倾角
                if(fabs(Ctrl->M.rake) > 180 ||
                   (Ctrl->M.full && (Ctrl->M.strike < 0 || Ctrl->M.strike > 360 || Ctrl->M.dip < 0 || Ctrl->M.dip > 90))) {
                    GRTRaiseError("Invalid receiver mechanism.");
                }
                break;
            }
            GRT_Common_Options_in_Switch((char)optopt);
        }
    }

    // 必须提供应力结果目录
    GRTCheckOptionActive(Ctrl, G);
    if(optind != argc) {
        GRTRaiseError("Unexpected positional argument %s. Use '-h' for help.", argv[optind]);
    }

    // 统一接收机制和逐点机制文件不能同时使用
    if(Ctrl->M.active && Ctrl->Q.active) {
        GRTRaiseError("-M and -Q are mutually exclusive.");
    }
}

/** 子模块主函数 */
int sproj_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    // 展开接收目录，按需读取替换接收机制的逐点文件
    size_t nr;
    char **dirs = grt_syn_output_receiver_directories(Ctrl->G.path, &nr);
    size_t npoints = 0;
    RCV_POINT *rcvs = Ctrl->Q.active ? grt_rcv_points_from_file(Ctrl->Q.path, &npoints) : NULL;
    if(rcvs && !rcvs[0].has_mechanism) {
        GRTRaiseError("-Q requires six columns.");
    }

    // 使用 -Q 时先核对所有接收位置，避免坐标不匹配时已经写出部分结果
    if(rcvs) {
        if(npoints != nr) {
            GRTRaiseError("Receiver count in -Q (%zu) does not match input directories (%zu).", npoints, nr);
        }

        // 按已确定的接收顺序逐点核对类型和坐标，有限接收断层不能使用 -Q 替换几何
        for(size_t ir = 0; ir < nr; ++ir) {
            char *path = NULL;
            GRT_SAFE_ASPRINTF(&path, "%s/stress_ZZ.sac", dirs[ir]);
            SACTRACE *trace = grt_read_SACTRACE(path, true);
            GRT_SAFE_FREE_PTR(path);
            const SACHEAD *hd = &trace->hd;
            if(GRT_SACHEAD_GET_RCV_FAULT_INDEX(hd) >= 0) {
                GRTRaiseError("Finite receivers require -Mrake[+f], without -Q.");
            }

            // 使用 SAC 头段中的坐标逐点比较，避免目录名两位小数的舍入影响
            if(!GRT_ISCLOSE(rcvs[ir].north, GRT_SACHEAD_GET_RCV_NORTH(hd)) ||
               !GRT_ISCLOSE(rcvs[ir].east, GRT_SACHEAD_GET_RCV_EAST(hd)) || !GRT_ISCLOSE(rcvs[ir].depth, -hd->stel * 1e-3)) {
                GRTRaiseError("Receiver coordinates in -Q at point %zu do not match %s.", ir, dirs[ir]);
            }
            grt_free_SACTRACE(trace);
        }
    }

    // 每个目录只读取一次应力，直接使用当前接收点的机制进行投影
    const char *stress_names[2][6] = {{"ZZ", "ZR", "ZT", "RR", "RT", "TT"}, {"ZZ", "ZN", "ZE", "NN", "NE", "EE"}};
    for(size_t ir = 0; ir < nr; ++ir) {
        char *path = NULL;
        GRT_SAFE_ASPRINTF(&path, "%s/stress_NN.sac", dirs[ir]);
        bool zne = access(path, F_OK) == 0;
        GRT_SAFE_FREE_PTR(path);

        // 根据标志性分量选择坐标系，再读取六个独立的应力分量
        SACTRACE *traces[6];
        for(int c = 0; c < 6; ++c) {
            GRT_SAFE_ASPRINTF(&path, "%s/stress_%s.sac", dirs[ir], stress_names[zne][c]);
            traces[c] = grt_read_SACTRACE(path, false);
            const SACHEAD *hd = &traces[c]->hd;
            if(hd->npts <= 0 || hd->delta <= 0 ||
               hd->npts != traces[0]->hd.npts || hd->delta != traces[0]->hd.delta || hd->b != traces[0]->hd.b) {
                GRTRaiseError("Invalid or inconsistent stress time sampling at %s.", path);
            }
            GRT_SAFE_FREE_PTR(path);
        }

        // SAC 头段作为默认机制，有限接收断层保留其走向和倾角，只允许调整滑动角
        const SACHEAD *hd = &traces[0]->hd;
        real_t strike = GRT_SACHEAD_GET_RCV_STRIKE(hd);
        real_t dip = GRT_SACHEAD_GET_RCV_DIP(hd);
        real_t rake = GRT_SACHEAD_GET_RCV_RAKE(hd);
        if(GRT_SACHEAD_GET_RCV_FAULT_INDEX(hd) >= 0) {
            if(Ctrl->M.active && Ctrl->M.full) {
                GRTRaiseError("Finite receivers require -Mrake[+f], without -Q.");
            }
            if(Ctrl->M.active && (Ctrl->M.force || fabs(rake) > 180)) {
                rake = Ctrl->M.rake;
            }
        } else if(rcvs) {
            // 普通接收点按文件中的接收顺序替换完整机制
            strike = rcvs[ir].strike;
            dip = rcvs[ir].dip;
            rake = rcvs[ir].rake;
        } else if(Ctrl->M.active) {
            // 普通接收点的手动机制必须完整，且不能带有限接收断层专用的 +f
            if(!Ctrl->M.full || Ctrl->M.force) {
                GRTRaiseError("Ordinary receivers require -Mstrike/dip/rake.");
            }
            strike = Ctrl->M.strike;
            dip = Ctrl->M.dip;
            rake = Ctrl->M.rake;
        }

        // 在生成输出前确认最终使用的机制完整且有效
        if(!rcvs && !Ctrl->M.full && (strike < 0 || strike > 360 || dip < 0 || dip > 90 || fabs(rake) > 180)) {
            GRTRaiseError("Undefined/invalid receiver mechanism at %s.", dirs[ir]);
        }

        SACTRACE *normal = grt_copy_SACTRACE(traces[0], true);
        SACTRACE *shear = grt_copy_SACTRACE(traces[0], true);

        // 在统一的全局坐标系中构造法向和滑动方向，逐时刻投影应力张量
        real_t nvec[3], tvec[3];
        grt_fault_plane_vectors(strike, dip, rake, nvec, tvec);
        for(int n = 0; n < normal->hd.npts; ++n) {
            real_t stress[6], sn, ts;
            for(int c = 0; c < 6; ++c) {
                stress[c] = traces[c]->data[n];
            }
            if(!zne) {
                grt_rot_zxy2zrt_symtensor2odr(-normal->hd.az * DEG1, stress);
            }
            grt_project_stress_to_fault_plane(stress, nvec, tvec, &sn, &ts);
            normal->data[n] = sn;
            shear->data[n] = ts;
        }

        // 投影结果继承接收信息，并记录实际采用的机制
        SACTRACE *outputs[2] = {normal, shear};
        const char *output_names[2] = {"sigma_n", "tau_s"};
        for(int c = 0; c < 2; ++c) {
            GRT_SACHEAD_SET_RCV_STRIKE(&outputs[c]->hd, strike);
            GRT_SACHEAD_SET_RCV_DIP(&outputs[c]->hd, dip);
            GRT_SACHEAD_SET_RCV_RAKE(&outputs[c]->hd, rake);
            snprintf(outputs[c]->hd.kcmpnm, sizeof(outputs[c]->hd.kcmpnm), "%s", output_names[c]);

            char *path = NULL;
            GRT_SAFE_ASPRINTF(&path, "%s/%s.sac", dirs[ir], output_names[c]);
            grt_write_SACTRACE(path, outputs[c]);
            GRT_SAFE_FREE_PTR(path);
            grt_free_SACTRACE(outputs[c]);
        }

        // 当前接收点计算完成，释放输入应力和目录路径
        for(int c = 0; c < 6; ++c) {
            grt_free_SACTRACE(traces[c]);
        }
        GRT_SAFE_FREE_PTR(dirs[ir]);
    }
    GRT_SAFE_FREE_PTR(dirs);
    GRT_SAFE_FREE_PTR(rcvs);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
