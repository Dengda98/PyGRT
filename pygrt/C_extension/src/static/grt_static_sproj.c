/**
 * @file   grt_static_sproj.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08-24
 *
 *    Project a static stress tensor onto receiver-fault geometry and write
 *    sigma_n and tau_s back to the input NetCDF file
 *
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    // 输入静态合成文件
    struct {
        bool active;
        char *s_ingrid;
    } G;
    // 接收断层形态
    struct {
        bool active;
        real_t strike;
        real_t dip;
        real_t rake;
        bool force_rake;
        bool has_geometry;
    } M;
    // 逐点接收断层形态文件
    struct {
        bool active;
        char *s_path;
    } Q;
} GRT_MODULE_CTRL;

/**
 * 释放参数控制结构体及其动态成员
 *
 * @param[in,out] Ctrl 参数控制结构体
 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    if(Ctrl == NULL) return;
    GRT_SAFE_FREE_PTR(Ctrl->G.s_ingrid);
    GRT_SAFE_FREE_PTR(Ctrl->Q.s_path);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印模块使用说明 */
static void print_help(void)
{
printf("\n"
"[grt static sproj] %s\n\n", GRT_VERSION);printf(
"    Project a static stress tensor onto receiver-fault geometry.\n"
"    The input file must contain the six stress_* variables produced by\n"
"    `static stress`. The result is written back to the same nc file as\n"
"    sigma_n and tau_s. Existing result variables are overwritten.\n"
"    Both grid and points layouts, and both ZNE and ZRT stress components,\n"
"    are supported.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt static sproj -G<ingrid> [-M<geometry>] [-Q<file>] [-h]\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<ingrid>    Input static synthesis file containing stress_* variables.\n"
"                  The result is written into this file.\n"
"\n"
"    -M<geometry>  For grid and ordinary points, use\n"
"                  <strike>/<dip>/<rake> in degree. For finite receiver\n"
"                  points, use <rake> or <rake>+f; +f forces the manual\n"
"                  rake for every point.\n"
"\n"
"    -Q<file>      For an ordinary points file, read replacement receiver\n"
"                  geometry. Each row must contain north east depth\n"
"                  strike dip rake, and rows must match the input point\n"
"                  order.\n"
"\n"
"    -h            Display this help message.\n"
"\n\n"
);
}


/**
 * 从命令行中读取模块选项
 *
 * @param[out] Ctrl 参数控制结构体
 * @param[in]  argc 命令行参数数量
 * @param[in]  argv 命令行参数数组
 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;

    // 解析输入文件、手动形态和逐点形态文件选项
    while((opt = getopt(argc, argv, ":G:M:Q:h")) != -1){
        switch(opt){
            case 'G':
                Ctrl->G.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->G.s_ingrid);
                Ctrl->G.s_ingrid = strdup(optarg);
                break;

            case 'M':
            {
                // 先按 + 拆分主参数和可选的强制标记
                char *string = strdup(optarg);
                char *token = strtok(string, "+");
                real_t strike = 0.0, dip = 0.0, rake = 0.0;
                int nscan = token == NULL ? 0 : sscanf(token, "%lf/%lf/%lf", &strike, &dip, &rake);
                if((nscan != 1) && (nscan != 3)){
                    GRT_SAFE_FREE_PTR(string);
                    GRTBadOptionError(M, "Expect <rake> or <strike>/<dip>/<rake>.");
                }

                Ctrl->M.active = true;
                Ctrl->M.has_geometry = nscan == 3;
                if(nscan == 3){
                    Ctrl->M.strike = strike;
                    Ctrl->M.dip = dip;
                    Ctrl->M.rake = rake;
                    if(!isfinite(strike) || (strike < 0.0) || (strike > 360.0)){
                        GRT_SAFE_FREE_PTR(string);
                        GRTBadOptionError(M, "Strike must be in [0, 360].");
                    }
                    if(!isfinite(dip) || (dip < 0.0) || (dip > 90.0)){
                        GRT_SAFE_FREE_PTR(string);
                        GRTBadOptionError(M, "Dip must be in [0, 90].");
                    }
                } else {
                    Ctrl->M.rake = rake;
                }
                if(!isfinite(rake) || (rake < -180.0) || (rake > 180.0)){
                    GRT_SAFE_FREE_PTR(string);
                    GRTBadOptionError(M, "Rake must be in [-180, 180].");
                }

                // 目前只支持 +f，用于强制所有有限接收点使用手动 rake
                Ctrl->M.force_rake = false;
                token = strtok(NULL, "+");
                if(token != NULL){
                    if(strcmp(token, "f") != 0 || (strtok(NULL, "+") != NULL)){
                        GRT_SAFE_FREE_PTR(string);
                        GRTBadOptionError(M, "Only the +f modifier is supported.");
                    }
                    Ctrl->M.force_rake = true;
                }
                GRT_SAFE_FREE_PTR(string);
                break;
            }

            case 'Q':
                Ctrl->Q.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->Q.s_path);
                Ctrl->Q.s_path = strdup(optarg);
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 输入文件必须通过 -G 指定，模块不接受位置参数
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, G);
    if(optind < argc){
        GRTRaiseError("Unexpected positional argument \"%s\". Use -G<ingrid>.\n", argv[optind]);
    }
}


/**
 * 规范化走向角到 [0, 360)
 *
 * @param[in] strike 原始走向角
 * @return 规范化后的走向角
 */
static real_t normalize_strike(real_t strike)
{
    strike = fmod(strike, 360.0);
    if(strike < 0.0) strike += 360.0;
    return strike;
}


/**
 * 检查一个接收断层几何是否定义且合法
 *
 * @param[in] strike 走向角
 * @param[in] dip    倾角
 * @param[in] rake   滑动角
 * @return 几何有效且完整时返回 true
 */
static bool geometry_is_defined(real_t strike, real_t dip, real_t rake)
{
    if(!isfinite(strike) || !isfinite(dip) || !isfinite(rake)) return false;
    if((dip < 0.0) || (dip > 90.0)) return false;
    if((rake == GRT_FINITE_FAULT_UNDEFINED_RAKE) ||
        (rake < -180.0) || (rake > 180.0)){
        return false;
    }
    return true;
}


/**
 * 检查从 NetCDF 读取的变量维度是否与接收布局一致
 *
 * @param[in] ncid            NetCDF 文件 ID
 * @param[in] varid           待检查的变量 ID
 * @param[in] name            变量名称
 * @param[in] expected_ndims  期望的维度数量
 * @param[in] expected_dimids 期望的维度 ID 数组
 */
static void check_var_dimensions(
    int ncid, int varid, const char *name, int expected_ndims,
    const int *expected_dimids)
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    if(ndims != expected_ndims){
        GRTRaiseError(
            "Variable \"%s\" has %d dimensions, expected %d.",
            name, ndims, expected_ndims);
    }
    NC_CHECK(nc_inq_vardimid(ncid, varid, dimids));
    for(int i = 0; i < expected_ndims; ++i){
        if(dimids[i] != expected_dimids[i]){
            GRTRaiseError("Variable \"%s\" has an unexpected dimension.", name);
        }
    }
}


/**
 * 查找可选变量，区分变量不存在与其他 NetCDF 错误
 *
 * @param[in]  ncid  NetCDF 文件 ID
 * @param[in]  name  变量名称
 * @param[out] varid 变量 ID
 * @return 找到变量时返回 true
 */
static bool find_optional_var(int ncid, const char *name, int *varid)
{
    int status = nc_inq_varid(ncid, name, varid);
    if(status == NC_NOERR) return true;
    if(status != NC_ENOTVAR) NC_CHECK(status);
    return false;
}


/**
 * 从属性或应力变量名判断输入应力的坐标形式
 *
 * @param[in] ncid NetCDF 文件 ID
 * @return 输入应力为 ZNE 坐标时返回 true
 */
static bool get_input_rot2ZNE(int ncid)
{
    int rot2ZNE;
    int status = nc_get_att_int(ncid, NC_GLOBAL, "rot2ZNE", &rot2ZNE);
    if(status == NC_NOERR) return rot2ZNE != 0;
    if(status != NC_ENOTATT) NC_CHECK(status);

    // 没有坐标属性时，根据应力变量名称判断坐标形式
    int varid;
    if(find_optional_var(ncid, "stress_ZN", &varid)) return true;
    if(find_optional_var(ncid, "stress_ZR", &varid)) return false;
    GRTRaiseError("Cannot determine whether input stress components are ZNE or ZRT.");
    return false;
}


/**
 * 读取六个独立应力分量
 *
 * @param[in]  ncid    NetCDF 文件 ID
 * @param[in]  channels 应力坐标分量名称
 * @param[in]  npts    接收点数量
 * @param[in]  ndims   应力变量维度数量
 * @param[in]  dimids  应力变量维度 ID 数组
 * @param[out] stress6 六个应力分量的点序列
 */
static void read_stress_components(
    int ncid, const char *channels, size_t npts,
    int ndims, const int *dimids, real_t *stress6)
{
    static const int component_pairs[6][2] = {
        {0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}
    };

    // 按固定顺序读取六个独立分量，便于后续逐点重组张量
    for(int i = 0; i < 6; ++i){
        char name[32];
        int varid;
        int c1 = component_pairs[i][0];
        int c2 = component_pairs[i][1];
        snprintf(name, sizeof(name), "stress_%c%c", channels[c1], channels[c2]);
        if(!find_optional_var(ncid, name, &varid)){
            GRTRaiseError("Input file does not contain variable \"%s\".", name);
        }
        check_var_dimensions(ncid, varid, name, ndims, dimids);
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, stress6 + i * npts));
    }
}


/**
 * 根据断层几何三要素获得法向矢量和滑动矢量
 *
 * @param[in]  strike 走向角
 * @param[in]  dip    倾角
 * @param[in]  rake   滑动角
 * @param[out] nvec   接收断层法向矢量
 * @param[out] tvec   rake 方向滑动矢量
 */
static void get_receiver_vectors(
    const real_t strike, const real_t dip, const real_t rake,
    real_t nvec[3], real_t tvec[3])
{
    // 角度转为弧度后，按 N、E、Z 顺序构造两个接收断层方向
    real_t stk = DEG1 * normalize_strike(strike);
    real_t dipp = DEG1 * dip;
    real_t rak = DEG1 * rake;

    real_t sdip = sin(dipp);
    real_t cdip = cos(dipp);
    real_t sstk = sin(stk);
    real_t cstk = cos(stk);
    real_t srak = sin(rak);
    real_t crak = cos(rak);

    // 矢量顺序为 N、E、Z
    nvec[0] = -sstk * sdip;
    nvec[1] = cstk * sdip;
    nvec[2] = cdip;

    tvec[0] = crak * cstk + srak * cdip * sstk;
    tvec[1] = crak * sstk - srak * cdip * cstk;
    tvec[2] = srak * sdip;
}


/**
 * 将一个 ZRT 应力张量转换为同一点的 ZNE 应力张量
 *
 * @param[in]     theta  R 轴相对 N 轴的方位角
 * @param[in,out] stress 应力张量分量
 */
static void convert_zrt_stress_to_zne(real_t theta, real_t stress[6])
{
    // theta 为 R 轴相对 N 轴的方位角，逆旋转得到 N/E 分量
    grt_rot_zxy2zrt_symtensor2odr(-theta, stress);
}


/**
 * 在 ZNE 坐标下将应力张量投影到接收断层的法向和滑动方向
 *
 * @param[in]  stress    ZNE 应力张量分量
 * @param[in]  nvec      接收断层法向矢量
 * @param[in]  tvec      rake 方向滑动矢量
 * @param[out] sigma_n   法向应力投影
 * @param[out] tau_s     rake 方向剪应力投影
 */
static void project_stress(
    const real_t stress[6], const real_t nvec[3], const real_t tvec[3],
    real_t *sigma_n, real_t *tau_s)
{
    // stress 顺序为 ZZ、ZN、ZE、NN、NE、EE，矢量顺序为 N、E、Z
    real_t traction[3];
    // 先计算法向量上的牵引力，再分别取法向和 rake 方向分量
    traction[0] = stress[3] * nvec[0] + stress[4] * nvec[1] + stress[1] * nvec[2];
    traction[1] = stress[4] * nvec[0] + stress[5] * nvec[1] + stress[2] * nvec[2];
    traction[2] = stress[1] * nvec[0] + stress[2] * nvec[1] + stress[0] * nvec[2];

    *sigma_n = traction[0] * nvec[0] + traction[1] * nvec[1] + traction[2] * nvec[2];
    // 剪应力沿接收断层滑动方向投影，正负由 rake 方向决定
    *tau_s = traction[0] * tvec[0] + traction[1] * tvec[1] + traction[2] * tvec[2];
}


/**
 * 读取普通 points 布局的 depth 变量
 *
 * @param[in] ncid       NetCDF 文件 ID
 * @param[in] recv_info  接收点布局信息
 * @return 接收点深度数组
 */
static real_t *read_point_depths(int ncid, const GRT_RECV_NC_INFO *recv_info)
{
    int depth_varid;
    // 深度只用于核对 -Q 文件中的点顺序，不参与应力投影
    real_t *depths = (real_t *)calloc(recv_info->npts, sizeof(real_t));
    if(!find_optional_var(ncid, "depth", &depth_varid)){
        GRT_SAFE_FREE_PTR(depths);
        GRTRaiseError("Points input file does not contain variable \"depth\".");
    }
    check_var_dimensions(ncid, depth_varid, "depth", 1, recv_info->dimids);
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, depth_varid, depths));
    return depths;
}


/**
 * 比较两个接收点坐标，允许文本写入 NetCDF 后的微小舍入误差
 *
 * @param[in] first  第一个坐标值
 * @param[in] second 第二个坐标值
 * @return 两个坐标值足够接近时返回 true
 */
static bool same_coordinate(real_t first, real_t second)
{
    real_t scale = fmax(1.0, fmax(fabs(first), fabs(second)));
    return fabs(first - second) <= 1e-8 * scale;
}


/**
 * 从 -Q 文件读取并检查逐点接收断层形态
 *
 * @param[in]  ncid       输入 NetCDF 文件 ID
 * @param[in]  recv_info  输入接收点布局信息
 * @param[in]  path       -Q 文件路径
 * @param[out] strikes    各接收点走向角
 * @param[out] dips       各接收点倾角
 * @param[out] rakes      各接收点滑动角
 */
static void load_geometry_from_Q(
    int ncid, const GRT_RECV_NC_INFO *recv_info, const char *path,
    real_t *strikes, real_t *dips, real_t *rakes)
{
    GRT_RECV_POINTS *q_points = grt_recv_points_from_file(path);
    if(!q_points->has_geometry){
        grt_recv_points_free(q_points);
        GRTRaiseError("-Q file \"%s\" must contain exactly 6 columns.", path);
    }
    if(q_points->npts != recv_info->npts){
        size_t q_npts = q_points->npts;
        grt_recv_points_free(q_points);
        GRTRaiseError(
            "-Q file \"%s\" has %zu points, but the input file has %zu points.",
            path, q_npts, recv_info->npts);
    }

    // 先读取输入文件深度，再逐点核对坐标、顺序和接收断层形态
    real_t *depths = read_point_depths(ncid, recv_info);
    for(size_t i = 0; i < recv_info->npts; ++i){
        if(!same_coordinate(q_points->norths[i], recv_info->norths[i]) ||
            !same_coordinate(q_points->easts[i], recv_info->easts[i]) ||
            !same_coordinate(q_points->depths[i], depths[i])){
            GRT_SAFE_FREE_PTR(depths);
            grt_recv_points_free(q_points);
            GRTRaiseError(
                "-Q file \"%s\" does not have the same point order and coordinates "
                "as the input file at point %zu.", path, i);
        }
        if(!geometry_is_defined(
            q_points->strikes[i], q_points->dips[i], q_points->rakes[i])){
            GRT_SAFE_FREE_PTR(depths);
            grt_recv_points_free(q_points);
            GRTRaiseError("Undefined or invalid receiver geometry in -Q at point %zu.", i);
        }
        strikes[i] = normalize_strike(q_points->strikes[i]);
        dips[i] = q_points->dips[i];
        rakes[i] = q_points->rakes[i];
    }
    GRT_SAFE_FREE_PTR(depths);
    grt_recv_points_free(q_points);
}


/**
 * 从普通 points 输入变量读取接收断层形态
 *
 * @param[in]  ncid       输入 NetCDF 文件 ID
 * @param[in]  recv_info  输入接收点布局信息
 * @param[in]  Ctrl       参数控制结构体
 * @param[out] strikes    各接收点走向角
 * @param[out] dips       各接收点倾角
 * @param[out] rakes      各接收点滑动角
 */
static void load_geometry_from_points(
    int ncid, const GRT_RECV_NC_INFO *recv_info,
    const GRT_MODULE_CTRL *Ctrl,
    real_t *strikes, real_t *dips, real_t *rakes)
{
    int strike_varid, dip_varid, rake_varid;
    // 普通 points 的三要素必须全部存在，并且使用 point 维
    bool has_strike = find_optional_var(ncid, "strike", &strike_varid);
    bool has_dip = find_optional_var(ncid, "dip", &dip_varid);
    bool has_rake = find_optional_var(ncid, "rake", &rake_varid);
    bool has_any_geometry = has_strike || has_dip || has_rake;
    bool has_all_geometry = has_strike && has_dip && has_rake;

    if(has_strike) check_var_dimensions(ncid, strike_varid, "strike", 1, recv_info->dimids);
    if(has_dip) check_var_dimensions(ncid, dip_varid, "dip", 1, recv_info->dimids);
    if(has_rake) check_var_dimensions(ncid, rake_varid, "rake", 1, recv_info->dimids);

    if(has_all_geometry){
        // 只有三要素完整时才从输入文件读取逐点形态
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, strike_varid, strikes));
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, dip_varid, dips));
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, rake_varid, rakes));
    }

    bool geometry_complete = has_all_geometry;
    if(has_all_geometry){
        for(size_t i = 0; i < recv_info->npts; ++i){
            if(!geometry_is_defined(strikes[i], dips[i], rakes[i])){
                geometry_complete = false;
                break;
            }
        }
    }

    if(Ctrl->M.active){
        // 手动形态优先于输入文件中的逐点形态
        if(has_any_geometry){
            GRTRaiseWarning(
                "The input point receiver geometry is ignored; using the new "
                "manual geometry from -M.");
        }
        for(size_t i = 0; i < recv_info->npts; ++i){
            strikes[i] = Ctrl->M.strike;
            dips[i] = Ctrl->M.dip;
            rakes[i] = Ctrl->M.rake;
        }
        return;
    }

    // 未设置 -M 时，必须使用文件中完整且有效的三要素
    if(!geometry_complete){
        GRTRaiseError(
            "Points input has no complete receiver geometry. Set -M<strike>/<dip>/<rake> "
            "or provide a 6-column -Q file.");
    }
    for(size_t i = 0; i < recv_info->npts; ++i){
        strikes[i] = normalize_strike(strikes[i]);
    }
}


/**
 * 从有限接收断层 points 输入变量读取接收断层形态
 *
 * @param[in]  ncid          输入 NetCDF 文件 ID
 * @param[in]  recv_info     输入接收点布局信息
 * @param[in]  nfault_dimid  nfault 维度 ID
 * @param[in]  Ctrl          参数控制结构体
 * @param[out] strikes       各接收点走向角
 * @param[out] dips          各接收点倾角
 * @param[out] rakes         各接收点滑动角
 */
static void load_geometry_from_finite_points(
    int ncid, const GRT_RECV_NC_INFO *recv_info,
    int nfault_dimid, const GRT_MODULE_CTRL *Ctrl,
    real_t *strikes, real_t *dips, real_t *rakes)
{
    size_t nfault;
    NC_CHECK(nc_inq_dimlen(ncid, nfault_dimid, &nfault));
    if(nfault == 0){
        GRTRaiseError("Finite receiver input has an empty nfault dimension.");
    }

    int strike_varid, dip_varid, rake_varid, offset_varid;
    if(!find_optional_var(ncid, "strike", &strike_varid) ||
        !find_optional_var(ncid, "dip", &dip_varid)){
        GRTRaiseError("Finite receiver input must contain strike and dip variables.");
    }
    if(!find_optional_var(ncid, "offset", &offset_varid)){
        GRTRaiseError("Finite receiver input must contain the offset variable.");
    }
    // 有限接收断层按 nfault 保存形态，再通过 offset 映射到 point
    check_var_dimensions(ncid, strike_varid, "strike", 1, &nfault_dimid);
    check_var_dimensions(ncid, dip_varid, "dip", 1, &nfault_dimid);
    check_var_dimensions(ncid, offset_varid, "offset", 1, &nfault_dimid);

    real_t *fault_strikes = (real_t *)calloc(nfault, sizeof(real_t));
    real_t *fault_dips = (real_t *)calloc(nfault, sizeof(real_t));
    real_t *fault_rakes = (real_t *)calloc(nfault, sizeof(real_t));
    int *offsets = (int *)calloc(nfault, sizeof(int));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, strike_varid, fault_strikes));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, dip_varid, fault_dips));
    NC_CHECK(nc_get_var_int(ncid, offset_varid, offsets));

    bool has_rake = find_optional_var(ncid, "rake", &rake_varid);
    if(has_rake){
        check_var_dimensions(ncid, rake_varid, "rake", 1, &nfault_dimid);
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, rake_varid, fault_rakes));
    }

    // 先统计未定义 rake，再判断是否必须使用或允许使用手动 rake
    bool any_undefined_rake = !has_rake;
    for(size_t ifault = 0; ifault < nfault; ++ifault){
        if(!isfinite(fault_strikes[ifault]) ||
            !isfinite(fault_dips[ifault]) ||
            (fault_dips[ifault] <= 0.0) || (fault_dips[ifault] > 90.0)){
            GRTRaiseError("Invalid strike or dip for finite receiver fault %zu.", ifault);
        }
        if(!has_rake) continue;
        if(fault_rakes[ifault] == GRT_FINITE_FAULT_UNDEFINED_RAKE ||
            !isfinite(fault_rakes[ifault])){
            any_undefined_rake = true;
        } else if((fault_rakes[ifault] < -180.0) || (fault_rakes[ifault] > 180.0)){
            GRTRaiseError("Invalid rake for finite receiver fault %zu.", ifault);
        }
    }

    if(!Ctrl->M.active && any_undefined_rake){
        GRT_SAFE_FREE_PTR(fault_strikes);
        GRT_SAFE_FREE_PTR(fault_dips);
        GRT_SAFE_FREE_PTR(fault_rakes);
        GRT_SAFE_FREE_PTR(offsets);
        GRTRaiseError(
            "Finite receiver input has undefined rake. Set -M<rake> to fill it, "
            "or use -M<rake>+f to override all rakes.");
    }
    if(Ctrl->M.active && !Ctrl->M.force_rake && !any_undefined_rake){
        GRT_SAFE_FREE_PTR(fault_strikes);
        GRT_SAFE_FREE_PTR(fault_dips);
        GRT_SAFE_FREE_PTR(fault_rakes);
        GRT_SAFE_FREE_PTR(offsets);
        GRTRaiseError(
            "All finite receiver rakes are already defined; use -M<rake>+f "
            "to force a manual rake.");
    }

    // 将每条断层的形态展开到其对应的连续 point 范围
    size_t start = 0;
    for(size_t ifault = 0; ifault < nfault; ++ifault){
        if((offsets[ifault] < 0) || ((size_t)offsets[ifault] <= start) ||
            ((size_t)offsets[ifault] > recv_info->npts)){
            GRTRaiseError("Invalid offset for finite receiver fault %zu.", ifault);
        }
        size_t end = (size_t)offsets[ifault];
        for(size_t i = start; i < end; ++i){
            strikes[i] = normalize_strike(fault_strikes[ifault]);
            dips[i] = fault_dips[ifault];
            if(Ctrl->M.active && Ctrl->M.force_rake){
                rakes[i] = Ctrl->M.rake;
            } else if(Ctrl->M.active &&
                (!has_rake ||
                 (fault_rakes[ifault] == GRT_FINITE_FAULT_UNDEFINED_RAKE) ||
                 !isfinite(fault_rakes[ifault]))){
                rakes[i] = Ctrl->M.rake;
            } else {
                rakes[i] = fault_rakes[ifault];
            }
        }
        start = end;
    }
    if(start != recv_info->npts){
        GRTRaiseError(
            "The last finite receiver offset (%zu) does not match point (%zu).",
            start, recv_info->npts);
    }

    GRT_SAFE_FREE_PTR(fault_strikes);
    GRT_SAFE_FREE_PTR(fault_dips);
    GRT_SAFE_FREE_PTR(fault_rakes);
    GRT_SAFE_FREE_PTR(offsets);
}


/**
 * 根据布局读取并确定每个接收点的断层形态
 *
 * @param[in]  ncid          输入 NetCDF 文件 ID
 * @param[in]  recv_info     输入接收点布局信息
 * @param[in]  Ctrl          参数控制结构体
 * @param[in]  finite_points 是否为有限接收断层 points 布局
 * @param[in]  nfault_dimid  nfault 维度 ID
 * @param[out] strikes       各接收点走向角
 * @param[out] dips          各接收点倾角
 * @param[out] rakes         各接收点滑动角
 */
static void load_receiver_geometry(
    int ncid, const GRT_RECV_NC_INFO *recv_info,
    const GRT_MODULE_CTRL *Ctrl, bool finite_points, int nfault_dimid,
    real_t *strikes, real_t *dips, real_t *rakes)
{
    if(finite_points){
        // nfault 存在时优先按有限接收断层规则处理
        if(Ctrl->Q.active){
            GRTRaiseError("-Q cannot be used with finite receiver points.");
        }
        if(Ctrl->M.active && Ctrl->M.has_geometry){
            GRTRaiseError("Finite receiver points require -M<rake> or -M<rake>+f.");
        }
        load_geometry_from_finite_points(
            ncid, recv_info, nfault_dimid, Ctrl,
            strikes, dips, rakes);
        return;
    }

    if(recv_info->layout == GRT_RECV_NC_LAYOUT_GRID){
        // grid 没有逐点形态，只能使用完整的手动形态
        if(Ctrl->Q.active){
            GRTRaiseError("-Q can only be used with a points-layout input file.");
        }
        if(!Ctrl->M.active || !Ctrl->M.has_geometry || Ctrl->M.force_rake){
            GRTRaiseError("Grid input has no receiver geometry; set -M<strike>/<dip>/<rake>.");
        }
        for(size_t i = 0; i < recv_info->npts; ++i){
            strikes[i] = Ctrl->M.strike;
            dips[i] = Ctrl->M.dip;
            rakes[i] = Ctrl->M.rake;
        }
        return;
    }

    if(Ctrl->Q.active && Ctrl->M.active){
        GRTRaiseError("Options -M and -Q are mutually exclusive for points input.");
    }
    if(Ctrl->Q.active){
        // -Q 提供普通 points 的逐点新形态
        load_geometry_from_Q(
            ncid, recv_info, Ctrl->Q.s_path, strikes, dips, rakes);
        return;
    }

    // 普通 points 在没有 -Q 时使用文件形态或统一的手动形态
    if(Ctrl->M.active && (!Ctrl->M.has_geometry || Ctrl->M.force_rake)){
        GRTRaiseError("Points input requires -M<strike>/<dip>/<rake>.");
    }
    load_geometry_from_points(
        ncid, recv_info, Ctrl, strikes, dips, rakes);
}


/**
 * 检查已有输出变量的维度
 *
 * @param[in]  ncid    NetCDF 文件 ID
 * @param[in]  name    输出变量名称
 * @param[in]  ndims   期望的维度数量
 * @param[in]  dimids  期望的维度 ID 数组
 * @param[out] varid   输出变量 ID
 * @return 变量存在且维度正确时返回 true
 */
static bool find_and_check_output_var(
    int ncid, const char *name, int ndims, const int *dimids, int *varid)
{
    if(!find_optional_var(ncid, name, varid)) return false;
    check_var_dimensions(ncid, *varid, name, ndims, dimids);
    return true;
}


/**
 * 将投影结果追加或覆盖写入原始 NetCDF 文件
 *
 * @param[in] ncid       输入 NetCDF 文件 ID
 * @param[in] ndims      输出变量维度数量
 * @param[in] dimids     输出变量维度 ID 数组
 * @param[in] sigma_n    法向应力投影
 * @param[in] tau_s       rake 方向剪应力投影
 */
static void write_projection_variables(
    int ncid, int ndims, const int *dimids,
    const real_t *sigma_n, const real_t *tau_s)
{
    int sigma_n_varid, tau_s_varid;
    // 已有结果变量必须与当前接收布局使用相同维度
    bool has_sigma_n = find_and_check_output_var(
        ncid, "sigma_n", ndims, dimids, &sigma_n_varid);
    bool has_tau_s = find_and_check_output_var(
        ncid, "tau_s", ndims, dimids, &tau_s_varid);

    if(has_sigma_n){
        GRTRaiseWarning("Variable \"sigma_n\" already exists and will be overwritten.");
    }
    if(has_tau_s){
        GRTRaiseWarning("Variable \"tau_s\" already exists and will be overwritten.");
    }

    if(!has_sigma_n || !has_tau_s){
        // 只定义缺失变量，避免改动已有变量的定义和数据
        NC_CHECK(nc_redef(ncid));
        if(!has_sigma_n){
            NC_CHECK(nc_def_var(ncid, "sigma_n", NC_REAL, ndims, dimids, &sigma_n_varid));
        }
        if(!has_tau_s){
            NC_CHECK(nc_def_var(ncid, "tau_s", NC_REAL, ndims, dimids, &tau_s_varid));
        }
        NC_CHECK(nc_enddef(ncid));
    }

    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, sigma_n_varid, sigma_n));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, tau_s_varid, tau_s));
}


/**
 * 执行静态应力投影模块
 *
 * @param[in] argc 命令行参数数量
 * @param[in] argv 命令行参数数组
 * @return 模块执行状态
 */
int static_sproj_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    GRTCheckFileExist(Ctrl->G.s_ingrid);

    // 以写模式打开输入文件，结果直接写回该文件
    int ncid;
    NC_CHECK(nc_open(Ctrl->G.s_ingrid, NC_WRITE, &ncid));

    GRT_RECV_NC_INFO recv_info;
    grt_recv_nc_info_load(ncid, &recv_info);
    int ndims = (recv_info.layout == GRT_RECV_NC_LAYOUT_POINTS) ? 1 : 2;
    size_t npts = recv_info.npts;
    // 通过 nfault 维度区分普通 points 和有限接收断层 points
    int nfault_dimid;
    int nfault_status = nc_inq_dimid(ncid, "nfault", &nfault_dimid);
    if((nfault_status != NC_NOERR) && (nfault_status != NC_EBADDIM)){
        NC_CHECK(nfault_status);
    }
    bool finite_points = (nfault_status == NC_NOERR);
    if((recv_info.layout == GRT_RECV_NC_LAYOUT_GRID) && finite_points){
        GRTRaiseError("Grid input must not contain an nfault dimension.");
    }

    bool rot2ZNE = get_input_rot2ZNE(ncid);
    const char *channels = rot2ZNE ? GRT_ZNE_CODES : GRT_ZRT_CODES;
    // 先把应力分量读入内存，避免在计算过程中反复访问 NetCDF
    real_t *stress6 = (real_t *)calloc(6 * npts, sizeof(real_t));
    read_stress_components(ncid, channels, npts, ndims, recv_info.dimids, stress6);

    real_t *strikes = (real_t *)calloc(npts, sizeof(real_t));
    real_t *dips = (real_t *)calloc(npts, sizeof(real_t));
    real_t *rakes = (real_t *)calloc(npts, sizeof(real_t));
    // 根据布局和命令行选项确定每个 point 的接收断层形态
    load_receiver_geometry(
        ncid, &recv_info, Ctrl, finite_points, nfault_dimid,
        strikes, dips, rakes);

    real_t *sigma_n = (real_t *)calloc(npts, sizeof(real_t));
    real_t *tau_s = (real_t *)calloc(npts, sizeof(real_t));
    // 逐点完成坐标转换、方向构造和应力投影
    for(size_t i = 0; i < npts; ++i){
        real_t stress[6] = {
            stress6[i],
            stress6[npts + i],
            stress6[2 * npts + i],
            stress6[3 * npts + i],
            stress6[4 * npts + i],
            stress6[5 * npts + i],
        };
        if(!rot2ZNE){
            // ZRT 应力需要根据当前接收点方位角转换为 ZNE
            real_t distance = hypot(recv_info.norths[i], recv_info.easts[i]);
            real_t theta = GRT_IS_ZERO(distance) ? 0.0 : atan2(recv_info.easts[i], recv_info.norths[i]);
            convert_zrt_stress_to_zne(theta, stress);
        }

        real_t nvec[3], tvec[3];
        get_receiver_vectors(strikes[i], dips[i], rakes[i], nvec, tvec);
        project_stress(stress, nvec, tvec, &sigma_n[i], &tau_s[i]);
    }

    // 追加缺失的结果变量，或覆盖已有结果变量
    write_projection_variables(ncid, ndims, recv_info.dimids, sigma_n, tau_s);
    NC_CHECK(nc_close(ncid));

    // 关闭文件后释放接收信息和计算缓存
    grt_recv_nc_info_free(&recv_info);
    GRT_SAFE_FREE_PTR(stress6);
    GRT_SAFE_FREE_PTR(strikes);
    GRT_SAFE_FREE_PTR(dips);
    GRT_SAFE_FREE_PTR(rakes);
    GRT_SAFE_FREE_PTR(sigma_n);
    GRT_SAFE_FREE_PTR(tau_s);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
