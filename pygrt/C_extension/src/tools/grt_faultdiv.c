/**
 * @file   grt_faultdiv.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-10
 *
 *    剖分 sCoulomb 格式有限断层，可设置时间函数与破裂延迟
 *
 */

#include "grt.h"

#undef I

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入断层 */
    struct {
        bool active;
        char *path;              ///< 输入 sCoulomb 文件路径
        size_t nfault;           ///< 原始断层数量
        FINITE_FAULT *faults;    ///< 原始断层及其剖分结果
    } F;
    /** 剖分尺寸 */
    struct {
        bool active;
        real_t dL;               ///< 沿走向目标剖分尺寸 (km)，各断层微调为整数等分
        real_t dW;               ///< 沿倾向目标剖分尺寸 (km)，各断层微调为整数等分
    } I;
    /** 破裂参数 */
    struct {
        bool active;
        bool has_o;              ///< 是否指定局部起裂点
        bool has_d;              ///< 是否指定传播方向
        real_t velocity;         ///< 恒定破裂速度 (km/s)
        real_t l;                ///< 起裂点沿走向坐标 (km)
        real_t w;                ///< 起裂点沿下倾坐标 (km)
        real_t degree;           ///< 类似 rake 的传播角度 (degree)
    } R;
    /** 时间函数 */
    struct {
        bool active;
        char *base;              ///< 不含整体延迟的时间函数选项值
        real_t delay;            ///< 整体时间偏移 (s)
    } D;
    /** 输出路径 */
    struct {
        bool active;
        char *path;              ///< 输出 sCoulomb 文件路径
    } O;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->F.path);
    grt_finite_fault_free(Ctrl->F.nfault, Ctrl->F.faults);
    GRT_SAFE_FREE_PTR(Ctrl->D.base);
    GRT_SAFE_FREE_PTR(Ctrl->O.path);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印使用说明 */
static void print_help()
{
printf("\n"
"[grt faultdiv] %s\n\n", GRT_VERSION);printf(
"    Subdivide sCoulomb Faults and Optionally Assign Rupture Delays\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt faultdiv -F<fault> -I<dL>/<dW> -O<outfile>\n"
"                 [-R<vel>[+o<along>/<down>][+d<degree>]]\n"
"                 [-D<tftype>[/<tfparams>][+d<delay>]] [-h]\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -F<fault>    Input sCoulomb-format fault file.\n"
"                 Only Kode=100/200/300 are supported. Preserve slips.\n"
"                 Ignore all input fields after \"bot\" column.\n"
"\n"
"    -I<dL>/<dW>  Positive target patch sizes (km): dL along strike,\n"
"                 dW down dip on the fault plane.\n"
"                 For each input rectangle, round length/dL and width/dW\n"
"                 to the nearest integers (at least 1, ties round up).\n"
"                 Adjust each size to span the full rectangle uniformly;\n"
"                 no smaller remainder patches at the boundaries.\n"
"\n"
"    -R<vel>[+o<along>/<down>][+d<degree>]\n"
"                 Positive rupture velocity (km/s). Choose one mode:\n"
"                 +o: radial rupture from local coordinates (km).\n"
"                     Require exactly one input fault record.\n"
"                     Origin is the start of the fault's top edge.\n"
"                     along: distance along strike toward the end.\n"
"                     down: distance down dip on the fault plane.\n"
"                 +d: planar front at a rake-like angle [-180,180],\n"
"                     0 along strike, +90 up dip, -90 down dip,\n"
"                     +/-180 opposite strike. Starts at the upstream\n"
"                     boundary of each complete input rectangle.\n"
"                     Multiple records rupture independently.\n"
"                 Rupture starts at t=0.\n"
"                 If -D is omitted, use an impulse (-Di).\n"
"\n"
"    -D<tftype>[/<tfparams>][+d<delay>]\n"
"                 Time function: i (impulse), p/T (parabola),\n"
"                 t/T1/T2/T3 (trapezoid), c/T1/T2 (cosine),\n"
"                 r/f0 (Ricker), 0/file (custom). Durations in s.\n"
"                 +d: nonnegative overall offset (s), added to\n"
"                     propagation time at each patch center.\n"
"                 Without -R, assign the same STF and offset to all.\n"
"                 Without -R and -D, write only 11 numeric columns.\n"
"                 Time function and custom paths are written as supplied.\n"
"\n"
"    -O<outfile>  Output subdivision fault file.\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt faultdiv -Ffault.inp -I1/1 -Ofault_subdiv.inp\n"
"    grt faultdiv -Ffault.inp -I1/1 -R2.5+o4/2 -Dc/1/1+d2 -Ofault_rupture.inp\n"
"    grt faultdiv -Ffault.inp -I1/1 -R2.5+d120 -Ofault_direction.inp\n"
"\n\n"
);
}


/** 从命令行中读取选项，处理后记录到参数控制结构体 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;
    char extra;
    while((opt = getopt(argc, argv, ":F:I:R:D:O:h")) != -1) {
        switch(opt) {
            // -F<fault>
            case 'F':
                Ctrl->F.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->F.path);
                Ctrl->F.path = strdup(optarg);
                break;

            // -I<dL>/<dW>
            case 'I':
                Ctrl->I.active = true;
                if(sscanf(optarg, "%lf/%lf%c", &Ctrl->I.dL, &Ctrl->I.dW, &extra) != 2 ||
                   Ctrl->I.dL <= 0.0 || Ctrl->I.dW <= 0.0) {
                    GRTBadOptionError(I, "Require two positive patch sizes.");
                }
                break;

            // -R<vel>[+o<along>/<down>][+d<degree>]
            case 'R': {
                Ctrl->R.active = true;
                Ctrl->R.has_o = Ctrl->R.has_d = false;
                char *cursor;
                Ctrl->R.velocity = strtod(optarg, &cursor);
                if(cursor == optarg || Ctrl->R.velocity <= 0.0) {
                    GRTBadOptionError(R, "Rupture velocity must be positive.");
                }
                while(*cursor) {
                    int nread = 0;
                    if(strncmp(cursor, "+o", 2) == 0 && !Ctrl->R.has_o) {
                        if(sscanf(cursor + 2, "%lf/%lf%n", &Ctrl->R.l, &Ctrl->R.w, &nread) != 2 ||
                           Ctrl->R.l < 0.0 || Ctrl->R.w < 0.0) {
                            GRTBadOptionError(R, "Require nonnegative local origin coordinates.");
                        }
                        Ctrl->R.has_o = true;
                    } else if(strncmp(cursor, "+d", 2) == 0 && !Ctrl->R.has_d) {
                        if(sscanf(cursor + 2, "%lf%n", &Ctrl->R.degree, &nread) != 1 || fabs(Ctrl->R.degree) > 180.0) {
                            GRTBadOptionError(R, "Direction must be in [-180, 180].");
                        }
                        Ctrl->R.has_d = true;
                    } else {
                        GRTBadOptionError(R, "Unsupported or repeated rupture modifier.");
                    }
                    cursor += 2 + nread;
                }
                if(Ctrl->R.has_o == Ctrl->R.has_d) {
                    GRTBadOptionError(R, "Specify exactly one of +o and +d.");
                }
                break;
            }

            // -D<tftype>[/<tfparams>][+d<delay>]
            case 'D': {
                Ctrl->D.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->D.base);
                Ctrl->D.base = strdup(optarg);
                Ctrl->D.delay = 0.0;
                // 按 +d 分出时间函数与整体延迟，保留科学计数法和自定义路径中的其他 + 字符
                char *suffix = strstr(Ctrl->D.base, "+d");
                if(suffix) {
                    *suffix = '\0';
                    if(sscanf(suffix + 1, "d%lf%c", &Ctrl->D.delay, &extra) != 1 || Ctrl->D.delay < 0.0) {
                        GRTBadOptionError(D, "Overall delay must be nonnegative.");
                    }
                }
                if(!*Ctrl->D.base) {
                    GRTBadOptionError(D, "Invalid time-function specification.");
                }
                break;
            }

            // -O<outfile>
            case 'O':
                Ctrl->O.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->O.path);
                Ctrl->O.path = strdup(optarg);
                break;

            GRT_Common_Options_in_Switch(optopt);
        }
    }

    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, F);
    GRTCheckOptionActive(Ctrl, I);
    GRTCheckOptionActive(Ctrl, O);
    if(Ctrl->R.active && !Ctrl->D.active) {
        Ctrl->D.base = strdup("i");
    }

    GRTCheckFileExist(Ctrl->F.path);
    if(grt_is_same_file(Ctrl->F.path, Ctrl->O.path)) {
        GRTBadOptionError(O, "Output must differ from the input fault file.");
    }
}


/** 子模块主函数 */
int faultdiv_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    // 以静态震源方式读入，只解析前 11 个数值列，忽略行末时间函数及其他字段
    Ctrl->F.faults = grt_finite_fault_load_coulomb(Ctrl->F.path, &Ctrl->F.nfault, true, 0.0, 1.0, NULL);
    if(Ctrl->R.has_o && Ctrl->F.nfault != 1) {
        GRTBadOptionError(R, "+o requires exactly one input fault record.");
    }

    for(size_t i = 0; i < Ctrl->F.nfault; ++i) {
        if(!KODE_IS_FINITE(Ctrl->F.faults[i].kode)) {
            GRTRaiseError("Only Kode=100/200/300 are supported (record %zu).", i + 1);
        }
    }

    // 先完成全部几何及起裂点检查，再打开输出文件
    for(size_t i = 0; i < Ctrl->F.nfault; ++i) {
        FINITE_FAULT *f = &Ctrl->F.faults[i];
        real_t L = hypot(f->east_end - f->east_begin, f->north_end - f->north_begin);
        real_t W = (f->bot - f->top) / sin(f->dip * DEG1);
        if(Ctrl->R.has_o && (Ctrl->R.l > L + 1e-10 || Ctrl->R.w > W + 1e-10)) {
            GRTBadOptionError(R, "Origin must lie inside the input fault.");
        }
        // 尺寸与目标间隔之比四舍五入，确定各方向的整数等分数，最少保留一块
        real_t nL = GRT_MAX(1.0, round(L / Ctrl->I.dL));
        real_t nW = GRT_MAX(1.0, round(W / Ctrl->I.dW));
        // 间隔向上取相邻浮点值，避免通用剖分函数的 ceil 因舍入误差多生成一块
        real_t dL = nextafter(L / nL, INFINITY);
        real_t dW = nextafter(W / nW, INFINITY);
        grt_finite_fault_subdiv(f, dL, dW, 0, NULL);
    }

    FILE *output = GRTCheckOpenFile(Ctrl->O.path, "w");
    // 输出统一使用分量格式的固定表头，第 7、8 列按 Kode 解释
    const char *headers[] = {"#", "X-start", "Y-start", "X-fin", "Y-fin", "Kode", "value1", "value2", "dip", "top", "bot"};
    for(int iheader = 0; iheader < 2; ++iheader) {
        for(size_t j = 0; j < sizeof(headers) / sizeof(*headers); ++j) {
            const char *label = iheader == 0 ? headers[j] : "---------------";
            int width = j == 0 ? 4 : (j == 5 ? 5 : 15);
            fprintf(output, "%s%*.*s", j ? " " : "", width, width, label);
        }
        fprintf(output, "\n");
    }

    real_t cosr = cos(Ctrl->R.degree * DEG1);  // 破裂方向角的余弦，沿正走向的分量
    real_t sinr = sin(Ctrl->R.degree * DEG1);  // 破裂方向角的正弦，沿上倾的分量
    const unsigned int ID = 1;  // 第一列为占位符，所有输出记录均填 1
    for(size_t i = 0; i < Ctrl->F.nfault; ++i) {
        FINITE_FAULT *f = &Ctrl->F.faults[i];
        real_t L = hypot(f->east_end - f->east_begin, f->north_end - f->north_begin);
        real_t sind = sin(f->dip * DEG1);
        real_t cosd = cos(f->dip * DEG1);
        real_t W = (f->bot - f->top) / sind;
        real_t alongE = (f->east_end - f->east_begin) / L;
        real_t alongN = (f->north_end - f->north_begin) / L;
        real_t upstream = GRT_MIN(0.0, L * cosr) + GRT_MIN(0.0, -W * sinr);
        // 使用已解析的滑动分量，rake/net slip 输入也统一转为分量格式
        real_t value1 = f->kode == KODE_TENSILE_REVERSE ? f->tensile : f->right_lateral;
        real_t value2 = f->kode == KODE_RTLAT_TENSILE ? f->tensile : f->reverse;
        for(size_t j = 0; j < f->nW * f->nL; ++j) {
            real_t halfL = 0.5 * f->length[j];
            real_t halfW = 0.5 * f->width[j];
            // 由子断层中心回推顶边，保留原始 Kode、分量及倾角
            real_t east0 = f->east[j] - halfL * alongE - halfW * cosd * alongN;
            real_t north0 = f->north[j] - halfL * alongN + halfW * cosd * alongE;
            real_t east1 = east0 + f->length[j] * alongE;
            real_t north1 = north0 + f->length[j] * alongN;
            // 中心回推的深度限制在母断层内，避免舍入误差将地表顶深写为负数
            real_t top = GRT_MAX(f->top, f->depth[j] - halfW * sind);
            real_t bot = GRT_MIN(f->bot, f->depth[j] + halfW * sind);
            fprintf(output, "%4u %15.6f %15.6f %15.6f %15.6f %5u %15.6f %15.6f %15.6f %15.6f %15.6f",
                ID, east0, north0, east1, north1, f->kode, value1, value2, f->dip, top, bot);

            if(Ctrl->R.active || Ctrl->D.active) {
                real_t delay = Ctrl->D.delay;
                if(Ctrl->R.active) {
                    // 破裂延迟按断层面局部中心计算，方向前沿以完整断层边界为基准
                    real_t w = (f->depth[j] - f->top) / sind;
                    real_t l = (f->east[j] - f->east_begin) * alongE + (f->north[j] - f->north_begin) * alongN;
                    real_t distance = Ctrl->R.has_o ? hypot(l - Ctrl->R.l, w - Ctrl->R.w) :
                                                     GRT_MAX(0.0, l * cosr - w * sinr - upstream);
                    delay += distance / Ctrl->R.velocity;
                }
                fprintf(output, "  -D%s+d%.6f", Ctrl->D.base, delay);
            }
            fprintf(output, "\n");
        }
    }
    if(ferror(output) || fclose(output) != 0) {
        GRTRaiseError("Failed to write %s.", Ctrl->O.path);
    }

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
