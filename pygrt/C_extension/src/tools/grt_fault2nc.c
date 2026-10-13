/**
 * @file   grt_fault2nc.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-10
 *
 *    将规则排布的 sCoulomb 剪切子断层整理为断层面 NetCDF 网格
 *
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    struct {
        bool active;
        char *path;              ///< 输入 sCoulomb 文件路径
        size_t nfault;           ///< 子断层数量
        FINITE_FAULT *faults;    ///< 子断层几何及滑动信息
    } F;
    struct {
        bool active;
        char *path;              ///< 输出 NetCDF 文件路径
    } O;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    GRT_SAFE_FREE_PTR(Ctrl->F.path);
    grt_finite_fault_free(Ctrl->F.nfault, Ctrl->F.faults);
    GRT_SAFE_FREE_PTR(Ctrl->O.path);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印使用说明 */
static void print_help()
{
printf("\n"
"[grt fault2nc] %s\n\n", GRT_VERSION);printf(
"    Convert a Regular sCoulomb Fault Plane to a NetCDF Grid\n"
"\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt fault2nc -F<fault> -O<outfile> [-h]\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -F<fault>    sCoulomb file with only Kode=100 rectangles.\n"
"                 Patches must form one complete plane, with equal\n"
"                 sizes and spacing along strike and down dip.\n"
"                 Input rows may be in any order.\n"
"                 Delay is taken from the final -D...+d value (s),\n"
"                 or zero if absent. Time functions are not sampled.\n"
"\n"
"    -O<outfile>  NetCDF output. Global attributes include strike/dip,\n"
"                 plane top/bot depths and patch dL/dW sizes (km).\n"
"                 Grid dimensions are width (down dip) and length\n"
"                 (along strike), with center coordinates in km.\n"
"\n"
"    -h           Display this help message.\n"
"\n"
"Example:\n"
"----------------------------------------------------------------\n"
"    grt fault2nc -Ffault_rupture.inp -Ofault.nc\n"
"\n");
}


/** 从命令行中读取选项，处理后记录到参数控制结构体 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;
    while((opt = getopt(argc, argv, ":F:O:h")) != -1) {
        switch(opt) {
            case 'F':
                Ctrl->F.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->F.path);
                Ctrl->F.path = strdup(optarg);
                break;

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
    GRTCheckOptionActive(Ctrl, O);

    GRTCheckFileExist(Ctrl->F.path);
    if(grt_is_same_file(Ctrl->F.path, Ctrl->O.path)) {
        GRTBadOptionError(O, "Output must differ from the input fault file.");
    }
}


/** 子模块主函数 */
int fault2nc_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    // 按震源读取几何、滑动分量和破裂延迟，dt<0 保留原始延迟且不采样时间函数
    Ctrl->F.faults = grt_finite_fault_load_coulomb(Ctrl->F.path, &Ctrl->F.nfault, true, -1.0, 1.0, NULL);
    size_t n = Ctrl->F.nfault;

    // 以第一块断层的顶边起点为参考原点，建立走向、倾向和法向坐标
    FINITE_FAULT *first = &Ctrl->F.faults[0];
    real_t strike = first->strike, dip = first->dip;
    real_t alongE = sin(strike * DEG1);
    real_t alongN = cos(strike * DEG1);
    real_t sind = sin(dip * DEG1);
    real_t cosd = cos(dip * DEG1);
    real_t dL = 0.0, dW = 0.0;
    real_t top = first->top;
    real_t bot = first->bot;

    // 每条记录仅取一个中心点，以平均长宽作为网格间距，减小文本舍入的影响
    for(size_t i = 0; i < n; ++i) {
        FINITE_FAULT *f = &Ctrl->F.faults[i];
        if(f->kode != KODE_RTLAT_REVERSE) {
            GRTRaiseError("Patches require Kode=100 (record %zu).", i + 1);
        }
        grt_finite_fault_subdiv(f, 0.0, 0.0, 0, NULL);
        dL += f->length[0] / n;
        dW += f->width[0] / n;
        top = GRT_MIN(top, f->top);
        bot = GRT_MAX(bot, f->bot);
    }
    // 允许文本坐标的约米级舍入误差，尺寸额外允许 0.01% 的相对误差
    real_t tol = 1e-3 + 1e-4 * GRT_MAX(dL, dW);
    real_t *localL = GRT_SAFE_CALLOC(n, sizeof(*localL));
    real_t *localW = GRT_SAFE_CALLOC(n, sizeof(*localW));
    real_t minL = INFINITY, minW = INFINITY;

    // 检查每块断层的几何，并将中心点投影到同一个断层面的局部坐标系
    for(size_t i = 0; i < n; ++i) {
        FINITE_FAULT *f = &Ctrl->F.faults[i];
        real_t east = f->east_begin - first->east_begin;
        real_t north = f->north_begin - first->north_begin;
        real_t depth = f->top - first->top;

        // 顶边起点到参考平面的有符号法向距离 (km)，参考单位法向量为 (-alongN*sind, alongE*sind, cosd)
        real_t normal = (-east * alongN + north * alongE) * sind + depth * cosd;

        // 当前断层顶边从起点指向终点的水平向量，edgeE、edgeN 分别为 east、north 分量 (km)
        real_t edgeE = f->east_end - f->east_begin;
        real_t edgeN = f->north_end - f->north_begin;

        // 沿顶边从起点到终点的法向距离增量 (km)
        real_t edge_normal = (-edgeE * alongN + edgeN * alongE) * sind;

        // 当前顶边单位向量与参考走向单位向量的点积，即两个走向夹角的余弦
        real_t orientation = (edgeE * alongE + edgeN * alongN) / f->length[0];

        // 沿当前断层下倾边从顶边到下边界的法向距离增量 (km)
        real_t bottom_normal = f->width[0] * (sin(f->dip * DEG1) * cosd - cos(f->dip * DEG1) * sind * orientation);

        // 四个顶点均须位于参考平面，顶边向量同时约束走向和长度
        // 四个法向距离依次对应顶边起点、顶边终点、下边界起点、下边界终点
        bool coplanar         = fabs(normal) <= tol && fabs(normal + edge_normal) <= tol &&
                                fabs(normal + bottom_normal) <= tol && fabs(normal + edge_normal + bottom_normal) <= tol;
        bool same_dip         = fabs(f->dip - dip) <= 0.05;
        bool same_strike_edge = hypot(edgeE - dL * alongE, edgeN - dL * alongN) <= tol;
        bool same_width       = fabs(f->width[0] - dW) <= tol;
        if(!coplanar || !same_dip || !same_strike_edge || !same_width) {
            GRTRaiseError("Patches must be coplanar with equal sizes and orientation (record %zu).", i + 1);
        }

        // 最小中心坐标用于定位网格，输出坐标再平移至完整断层的上角点
        east = f->east[0] - first->east_begin;
        north = f->north[0] - first->north_begin;
        depth = f->depth[0] - first->top;
        localL[i] = east * alongE + north * alongN;
        localW[i] = (east * alongN - north * alongE) * cosd + depth * sind;
        minL = GRT_MIN(minL, localL[i]);
        minW = GRT_MIN(minW, localW[i]);
    }

    // 按局部坐标定位网格单元，所有单元必须存在且仅出现一次
    size_t nL = 0, nW = 0;
    for(size_t i = 0; i < n; ++i) {
        real_t iL = round((localL[i] - minL) / dL);
        real_t iW = round((localW[i] - minW) / dW);
        bool regular_strike = fabs(localL[i] - minL - iL * dL) <= tol;
        bool regular_dip    = fabs(localW[i] - minW - iW * dW) <= tol;
        bool valid_index    = iL < n && iW < n;
        if(!regular_strike || !regular_dip || !valid_index) {
            GRTRaiseError("Patch centers must be regularly spaced (record %zu).", i + 1);
        }
        nL = GRT_MAX(nL, (size_t)iL + 1);
        nW = GRT_MAX(nW, (size_t)iW + 1);
    }
    // 先检查矩形单元总数，再检查重复记录，二者共同保证没有缺块
    if(nL > n / nW || nL * nW != n) {
        GRTRaiseError("Patches must cover a complete rectangular grid without gaps.");
    }
    size_t *order = GRT_SAFE_CALLOC(n, sizeof(*order));  // 保存原始记录索引加一，零表示单元尚未填充
    for(size_t i = 0; i < n; ++i) {
        size_t iL = (size_t)round((localL[i] - minL) / dL);
        size_t iW = (size_t)round((localW[i] - minW) / dW);
        size_t cell = iW * nL + iL;
        if(order[cell]) {
            GRTRaiseError("Duplicate patch at grid cell %zu/%zu.", iL, iW);
        }
        order[cell] = i + 1;
    }

    // 全部校验通过后创建文件，二维变量按 (width, length) 保存，沿走向变化最快
    int ncid, dims[2], coords[2], vars[6];
    const char *dimnames[2] = {"width", "length"};
    const char *varnames[6] = {"north", "east", "depth", "rake", "slip", "delay"};
    const char *units[6] = {"km", "km", "km", "degree", "m", "s"};
    size_t sizes[2] = {nW, nL};
    NC_CHECK(nc_create(Ctrl->O.path, NC_CLOBBER, &ncid));
    for(int k = 0; k < 2; ++k) {
        NC_CHECK(nc_def_dim(ncid, dimnames[k], sizes[k], &dims[k]));
        NC_CHECK(nc_def_var(ncid, dimnames[k], NC_REAL, 1, &dims[k], &coords[k]));
        NC_CHECK(nc_put_att_text(ncid, coords[k], "units", 2, "km"));
    }
    for(int k = 0; k < 6; ++k) {
        NC_CHECK(nc_def_var(ncid, varnames[k], NC_REAL, 2, dims, &vars[k]));
        NC_CHECK(nc_put_att_text(ncid, vars[k], "units", strlen(units[k]), units[k]));
    }
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "strike", NC_REAL, 1, &strike));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "dip", NC_REAL, 1, &dip));
    // 记录完整断层面的深度边界和平均子断层尺寸，便于直接设置绘图坐标轴
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "top", NC_REAL, 1, &top));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "bot", NC_REAL, 1, &bot));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "dL", NC_REAL, 1, &dL));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(ncid, NC_GLOBAL, "dW", NC_REAL, 1, &dW));
    NC_CHECK(nc_enddef(ncid));

    // 一维坐标从半个单元边长起算，表示相对于完整断层上角点的中心位置
    real_t *values = GRT_SAFE_CALLOC(n, sizeof(*values));
    for(int k = 0; k < 2; ++k) {
        for(size_t i = 0; i < sizes[k]; ++i) {
            values[i] = (i + 0.5) * (k == 0 ? dW : dL);
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, coords[k], values));
    }

    // 按网格顺序重排记录，从剪切源项读取实际滑动角，将单位面积矩势 (cm^3/km^2) 换回滑动量 (m)
    for(int k = 0; k < 6; ++k) {
        for(size_t cell = 0; cell < n; ++cell) {
            size_t i = order[cell] - 1;
            FINITE_FAULT *f = &Ctrl->F.faults[i];
            switch(k) {
                case 0: values[cell] = f->north[0]; break;
                case 1: values[cell] = f->east[0]; break;
                case 2: values[cell] = f->depth[0]; break;
                case 3: values[cell] = f->terms[0].mechanism[2]; break;
                case 4: values[cell] = f->terms[0].scale * 1e-12; break;
                case 5: values[cell] = f->stf_delay; break;
            }
        }
        NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, vars[k], values));
    }
    NC_CHECK(nc_close(ncid));
    GRT_SAFE_FREE_PTR(values);
    GRT_SAFE_FREE_PTR(order);
    GRT_SAFE_FREE_PTR(localL);
    GRT_SAFE_FREE_PTR(localW);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
