/**
 * @file   grt_xy2geo_geo2xy.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08-26
 *
 *    Convert between local north/east and geographic latitude/longitude
 *    coordinates for static NetCDF files and coordinate text files
 *
 */

#include "grt.h"

#include <errno.h>

#define EARTH_RADIUS_KM 6371.0

/** 该子模块的参数控制结构体 */
typedef struct {
    // 输入 NetCDF 文件
    struct {
        bool active;
        char *s_path;
    } G;
    // 输入文本坐标文件
    struct {
        bool active;
        char *s_path;
    } Q;
    // 输出文件
    struct {
        bool active;
        char *s_path;
    } O;
    // 参考点经纬度
    struct {
        bool active;
        real_t lat0;
        real_t lon0;
    } C;
} GRT_MODULE_CTRL;

/** 输入 NetCDF 文件中的坐标布局信息 */
typedef struct {
    GRT_RECV_NC_LAYOUT layout;
    size_t first_count;
    size_t second_count;
} GRT_COORD_INFO;

/** 坐标转换方向 */
typedef enum {
    GRT_XY2GEO = 0,
    GRT_GEO2XY
} GRT_TRANSFORM_DIRECTION;

/**
 * 释放参数控制结构体及其动态成员
 *
 * @param[in,out] Ctrl 参数控制结构体
 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    if(Ctrl == NULL) return;
    GRT_SAFE_FREE_PTR(Ctrl->G.s_path);
    GRT_SAFE_FREE_PTR(Ctrl->Q.s_path);
    GRT_SAFE_FREE_PTR(Ctrl->O.s_path);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/**
 * 打印模块使用说明
 *
 * @param[in] direction 坐标转换方向
 */
static void print_help_common(GRT_TRANSFORM_DIRECTION direction)
{
const bool to_geo = direction == GRT_XY2GEO;
const char *module = to_geo ? "xy2geo" : "geo2xy";
const char *description = to_geo
? "Convert local north/east coordinates to latitude/longitude."
: "Convert latitude/longitude coordinates to local north/east.";
const char *input_description = to_geo
? "local north/east coordinates"
: "geographic latitude/longitude coordinates";

printf("\n"
"[grt %s] %s\n\n", module, GRT_VERSION);
printf(
"    %s\n"
"    The input coordinates are in km for north/east and in degree\n"
"    for latitude/longitude. The reference is latitude/longitude\n"
"    in degree. Longitude is normalized to [-180, 180).\n"
"    Grid, ordinary points, finite receiver-point layouts, and\n"
"    coordinate text files are supported.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt %s -G<ingrid> -O<outgrid> -C<lat0>/<lon0> [-h]\n"
"    grt %s -Q<infile> -O<outfile> -C<lat0>/<lon0> [-h]\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -G<ingrid>    Input static NetCDF file containing %s.\n"
"                  The input file is not modified.\n"
"\n"
"    -Q<infile>    Input text file containing %s in the first two\n"
"                  columns. Empty lines and lines beginning with #\n"
"                  are copied without coordinate conversion.\n"
"\n"
"    -O<outfile>   Output NetCDF or text file path. The output file\n"
"                  must be different from the input file.\n"
"\n"
"    -C<lat0>/<lon0>\n"
"                  Reference latitude/longitude in degree. The\n"
"                  latitude must be in (-90, 90).\n"
"                  longitude must be in [-180, 180].\n"
"\n"
"    -h            Display this help message.\n"
"\n\n",
description, module, module, input_description, input_description);
}


/**
 * 读取参考点经纬度
 *
 * @param[out] Ctrl 参数控制结构体
 * @param[in]  arg  纬度/经度字符串
 */
static void parse_reference(GRT_MODULE_CTRL *Ctrl, const char *arg)
{
    real_t lat0, lon0;
    char extra;
    if(sscanf(arg, "%lf/%lf %c", &lat0, &lon0, &extra) != 2){
        GRTBadOptionError(C, "Expect <latitude>/<longitude>.");
    }
    if((!isfinite(lat0)) || (lat0 <= -90.0) || (lat0 >= 90.0)){
        GRTBadOptionError(C, "Latitude must be in (-90, 90).");
    }
    if((!isfinite(lon0)) || (lon0 < -180.0) || (lon0 > 180.0)){
        GRTBadOptionError(C, "Longitude must be in [-180, 180].");
    }
    Ctrl->C.active = true;
    Ctrl->C.lon0 = lon0;
    Ctrl->C.lat0 = lat0;
}


/**
 * 从命令行中读取模块选项
 *
 * @param[out] Ctrl      参数控制结构体
 * @param[in]  argc      命令行参数数量
 * @param[in]  argv      命令行参数数组
 * @param[in]  direction 坐标转换方向
 */
static void getopt_from_command(
    GRT_MODULE_CTRL *Ctrl, int argc, char **argv,
    GRT_TRANSFORM_DIRECTION direction)
{
    int opt;

    optind = 1;
    while((opt = getopt(argc, argv, ":G:Q:O:C:h")) != -1){
        switch(opt){
            case 'G':
                Ctrl->G.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->G.s_path);
                Ctrl->G.s_path = strdup(optarg);
                break;

            case 'Q':
                Ctrl->Q.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->Q.s_path);
                Ctrl->Q.s_path = strdup(optarg);
                break;

            case 'O':
                Ctrl->O.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->O.s_path);
                Ctrl->O.s_path = strdup(optarg);
                break;

            case 'C':
                parse_reference(Ctrl, optarg);
                break;

            case 'h':
                print_help_common(direction);
                exit(EXIT_SUCCESS);
                break;

            case ':':
                GRTMissArgsError((char)optopt, "");
                break;

            case '?':
            default:
                GRTInvalidOptionError((char)optopt, "");
                break;
        }
    }

    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, O);
    GRTCheckOptionActive(Ctrl, C);
    if((Ctrl->G.active && Ctrl->Q.active)){
        GRTRaiseError("Options \"-G\" and \"-Q\" are mutually exclusive.\n");
    }
    if((!Ctrl->G.active) && (!Ctrl->Q.active)){
        GRTRaiseError("Exactly one of options \"-G\" and \"-Q\" must be set.\n");
    }
    if(optind < argc){
        GRTRaiseError(
            "Unexpected positional argument \"%s\". Use -G or -Q, -O and -C.\n",
            argv[optind]);
    }
}


/**
 * 检查输出路径不是输入文件本身
 *
 * @param[in] input_path  输入文件路径
 * @param[in] output_path 输出文件路径
 */
static void check_output_path(const char *input_path, const char *output_path)
{
    struct stat input_stat, output_stat;
    if(strcmp(input_path, output_path) == 0){
        GRTRaiseError("Input and output files must be different.\n");
    }
    if(((stat(input_path, &input_stat) == 0) &&
        (stat(output_path, &output_stat) == 0)) &&
        ((input_stat.st_dev == output_stat.st_dev) &&
        (input_stat.st_ino == output_stat.st_ino))){
        GRTRaiseError("Input and output files must be different.\n");
    }
}


/**
 * 获取当前转换方向对应的坐标名称
 *
 * @param[in]  direction 坐标转换方向
 * @param[out] first     第一坐标名称
 * @param[out] second    第二坐标名称
 */
static void get_coordinate_names(
    GRT_TRANSFORM_DIRECTION direction, const char **first, const char **second)
{
    if(direction == GRT_XY2GEO){
        *first = "north";
        *second = "east";
    } else {
        *first = "lat";
        *second = "lon";
    }
}


/**
 * 检查坐标变量是一维变量，且使用指定的坐标维度
 *
 * @param[in] ncid           NetCDF 文件 ID
 * @param[in] varid          坐标变量 ID
 * @param[in] name           坐标变量名称
 * @param[in] expected_dimid 期望的坐标维度 ID
 */
static void check_coordinate_var(
    int ncid, int varid, const char *name, int expected_dimid)
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    NC_CHECK(nc_inq_varndims(ncid, varid, &ndims));
    if(ndims != 1){
        GRTRaiseError("Coordinate variable \"%s\" must be one-dimensional.", name);
    }
    NC_CHECK(nc_inq_vardimid(ncid, varid, dimids));
    if(dimids[0] != expected_dimid){
        GRTRaiseError("Coordinate variable \"%s\" has an unexpected dimension.", name);
    }
}


/**
 * 读取并检查输入文件的静态接收坐标布局
 *
 * @param[in]  ncid      文件 ID
 * @param[in]  direction 坐标转换方向
 * @param[out] info      坐标布局信息
 */
static void inspect_coordinate_layout(
    int ncid, GRT_TRANSFORM_DIRECTION direction, GRT_COORD_INFO *info)
{
    const char *first_name, *second_name;
    get_coordinate_names(direction, &first_name, &second_name);

    int first_varid, second_varid;
    NC_CHECK(nc_inq_varid(ncid, first_name, &first_varid));
    NC_CHECK(nc_inq_varid(ncid, second_name, &second_varid));

    info->layout = grt_recv_nc_get_layout(ncid);
    if(info->layout == GRT_RECV_NC_LAYOUT_GRID){
        int first_dimid, second_dimid;
        NC_CHECK(nc_inq_dimid(ncid, first_name, &first_dimid));
        NC_CHECK(nc_inq_dimid(ncid, second_name, &second_dimid));
        NC_CHECK(nc_inq_dimlen(ncid, first_dimid, &info->first_count));
        NC_CHECK(nc_inq_dimlen(ncid, second_dimid, &info->second_count));
        check_coordinate_var(ncid, first_varid, first_name, first_dimid);
        check_coordinate_var(ncid, second_varid, second_name, second_dimid);
        return;
    }

    int point_dimid;
    NC_CHECK(nc_inq_dimid(ncid, "point", &point_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, point_dimid, &info->first_count));
    info->second_count = info->first_count;
    check_coordinate_var(ncid, first_varid, first_name, point_dimid);
    check_coordinate_var(ncid, second_varid, second_name, point_dimid);
}


/**
 * 将经度规范化到 [-180, 180)
 *
 * @param[in] longitude 待规范化的经度
 * @return 规范化后的经度
 */
static real_t normalize_longitude(real_t longitude)
{
    longitude = fmod(longitude + 180.0, 360.0);
    if(longitude < 0.0) longitude += 360.0;
    return longitude - 180.0;
}


/**
 * 对一个坐标点执行局部切平面转换
 *
 * @param[in]  first     第一坐标
 * @param[in]  second    第二坐标
 * @param[in]  lon0      参考点经度
 * @param[in]  lat0      参考点纬度
 * @param[in]  direction 坐标转换方向
 * @param[out] out_first 转换后的第一坐标
 * @param[out] out_second 转换后的第二坐标
 */
static void transform_coordinate_pair(
    real_t first, real_t second, real_t lon0, real_t lat0,
    GRT_TRANSFORM_DIRECTION direction, real_t *out_first, real_t *out_second)
{
    const real_t km_per_lat_deg = EARTH_RADIUS_KM * DEG1;
    const real_t km_per_lon_deg = EARTH_RADIUS_KM * cos(lat0 * DEG1) * DEG1;

    if((!isfinite(first)) || (!isfinite(second))){
        GRTRaiseError("Coordinate values must be finite.");
    }

    if(direction == GRT_XY2GEO){
        real_t latitude = lat0 + first / km_per_lat_deg;
        real_t longitude = lon0 + second / km_per_lon_deg;
        if((!isfinite(latitude)) || (latitude < -90.0) || (latitude > 90.0)){
            GRTRaiseError("Converted latitude is outside [-90, 90].");
        }
        if((!isfinite(longitude))){
            GRTRaiseError("Converted longitude is not finite.");
        }
        *out_first = latitude;
        *out_second = normalize_longitude(longitude);
    } else {
        if((first < -90.0) || (first > 90.0)){
            GRTRaiseError("Latitude must be in [-90, 90].");
        }
        *out_first = (first - lat0) * km_per_lat_deg;
        *out_second = normalize_longitude(second - lon0) * km_per_lon_deg;
    }
}


/**
 * 转换输出文件中的坐标变量
 *
 * @param[in]  ncid      输出文件 ID
 * @param[in]  info      坐标布局信息
 * @param[in]  Ctrl      参数控制结构体
 * @param[in]  direction 坐标转换方向
 */
static void convert_output_coordinates(
    int ncid, const GRT_COORD_INFO *info,
    const GRT_MODULE_CTRL *Ctrl, GRT_TRANSFORM_DIRECTION direction)
{
    const char *first_name, *second_name;
    get_coordinate_names(direction, &first_name, &second_name);

    int first_varid, second_varid;
    NC_CHECK(nc_inq_varid(ncid, first_name, &first_varid));
    NC_CHECK(nc_inq_varid(ncid, second_name, &second_varid));

    real_t *firsts = (real_t *)calloc(info->first_count, sizeof(real_t));
    real_t *seconds = (real_t *)calloc(info->second_count, sizeof(real_t));
    if(((firsts == NULL) && (info->first_count != 0)) ||
        ((seconds == NULL) && (info->second_count != 0))){
        GRT_SAFE_FREE_PTR(firsts);
        GRT_SAFE_FREE_PTR(seconds);
        GRTRaiseError("Unable to allocate memory for coordinate variables.");
    }
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, first_varid, firsts));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, second_varid, seconds));

    for(size_t i = 0; i < info->first_count; ++i){
        real_t first, second;
        transform_coordinate_pair(
            firsts[i], 0.0, Ctrl->C.lon0, Ctrl->C.lat0,
            direction, &first, &second);
        firsts[i] = first;
    }
    for(size_t i = 0; i < info->second_count; ++i){
        real_t first, second;
        transform_coordinate_pair(
            direction == GRT_XY2GEO ? 0.0 : Ctrl->C.lat0, seconds[i],
            Ctrl->C.lon0, Ctrl->C.lat0, direction, &first, &second);
        seconds[i] = second;
    }

    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, first_varid, firsts));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, second_varid, seconds));

    GRT_SAFE_FREE_PTR(firsts);
    GRT_SAFE_FREE_PTR(seconds);
}


/**
 * 检查输出文件中是否已经存在目标名称
 *
 * @param[in] ncid  NetCDF 文件 ID
 * @param[in] name  目标名称
 * @param[in] is_dim 是否检查维度名称
 */
static void check_target_name_available(int ncid, const char *name, bool is_dim)
{
    int id;
    int status = is_dim ? nc_inq_dimid(ncid, name, &id) : nc_inq_varid(ncid, name, &id);
    int not_found = is_dim ? NC_EBADDIM : NC_ENOTVAR;
    if(status == NC_NOERR){
        GRTRaiseError("Output already contains the target name \"%s\".", name);
    }
    if(status != not_found){
        int error_status = status;
        NC_CHECK(error_status);
    }
}


/**
 * 写入参考点属性并重命名输出文件中的坐标维度和坐标变量
 *
 * @param[in]  ncid      输出文件 ID
 * @param[in]  layout    接收坐标布局
 * @param[in]  Ctrl      参数控制结构体
 * @param[in]  direction 坐标转换方向
 */
static void write_reference_and_rename_coordinates(
    int ncid, GRT_RECV_NC_LAYOUT layout, const GRT_MODULE_CTRL *Ctrl,
    GRT_TRANSFORM_DIRECTION direction)
{
    const char *first_name, *second_name;
    const char *target_first, *target_second;
    get_coordinate_names(direction, &first_name, &second_name);
    if(direction == GRT_XY2GEO){
        target_first = "lat";
        target_second = "lon";
    } else {
        target_first = "north";
        target_second = "east";
    }

    int first_varid, second_varid;
    NC_CHECK(nc_inq_varid(ncid, first_name, &first_varid));
    NC_CHECK(nc_inq_varid(ncid, second_name, &second_varid));
    check_target_name_available(ncid, target_first, false);
    check_target_name_available(ncid, target_second, false);

    int first_dimid = -1, second_dimid = -1;
    if(layout == GRT_RECV_NC_LAYOUT_GRID){
        NC_CHECK(nc_inq_dimid(ncid, first_name, &first_dimid));
        NC_CHECK(nc_inq_dimid(ncid, second_name, &second_dimid));
        check_target_name_available(ncid, target_first, true);
        check_target_name_available(ncid, target_second, true);
    }

    NC_CHECK(nc_redef(ncid));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(
        ncid, NC_GLOBAL, "lat0", NC_REAL, 1, &Ctrl->C.lat0));
    NC_CHECK(NC_FUNC_REAL(nc_put_att)(
        ncid, NC_GLOBAL, "lon0", NC_REAL, 1, &Ctrl->C.lon0));
    if(layout == GRT_RECV_NC_LAYOUT_GRID){
        NC_CHECK(nc_rename_dim(ncid, first_dimid, target_first));
        NC_CHECK(nc_rename_dim(ncid, second_dimid, target_second));
    }
    NC_CHECK(nc_rename_var(ncid, first_varid, target_first));
    NC_CHECK(nc_rename_var(ncid, second_varid, target_second));
    NC_CHECK(nc_enddef(ncid));
}


/**
 * 解析文本坐标行中的前两列
 *
 * @param[in]  line         文本行
 * @param[in]  lineno       行号
 * @param[in]  path         文件路径
 * @param[out] first        第一列坐标
 * @param[out] second       第二列坐标
 * @param[out] first_start  第一列数值起始位置
 * @param[out] first_end    第一列数值结束位置
 * @param[out] second_start 第二列数值起始位置
 * @param[out] second_end   第二列数值结束位置
 */
static void parse_text_coordinate_line(
    const char *line, size_t lineno, const char *path,
    real_t *first, real_t *second,
    const char **first_start, const char **first_end,
    const char **second_start, const char **second_end)
{
    const char *cursor = line;
    while(isspace((unsigned char)*cursor)) ++cursor;
    *first_start = cursor;

    errno = 0;
    char *end = NULL;
    *first = strtod(cursor, &end);
    if((end == cursor) || (errno == ERANGE) || (!isfinite(*first))){
        GRTRaiseError(
            "Invalid coordinate line %zu in \"%s\" "
            "(expect at least two finite numeric columns).", lineno, path);
    }
    *first_end = end;

    cursor = end;
    while(isspace((unsigned char)*cursor)) ++cursor;
    *second_start = cursor;
    errno = 0;
    *second = strtod(cursor, &end);
    if((end == cursor) || (errno == ERANGE) || (!isfinite(*second))){
        GRTRaiseError(
            "Invalid coordinate line %zu in \"%s\" "
            "(expect at least two finite numeric columns).", lineno, path);
    }
    *second_end = end;
}


/**
 * 转换文本坐标文件
 *
 * @param[in] input_path 输入文本文件路径
 * @param[in] output_path 输出文本文件路径
 * @param[in] Ctrl 参数控制结构体
 * @param[in] direction 坐标转换方向
 */
static void convert_text_file(
    const char *input_path, const char *output_path,
    const GRT_MODULE_CTRL *Ctrl, GRT_TRANSFORM_DIRECTION direction)
{
    FILE *in = fopen(input_path, "r");
    if(in == NULL){
        GRTFileOpenError(input_path);
    }
    FILE *out = fopen(output_path, "w");
    if(out == NULL){
        fclose(in);
        GRTFileOpenError(output_path);
    }

    if(fprintf(out, "# %.15g %.15g\n", Ctrl->C.lat0, Ctrl->C.lon0) < 0){
        GRTRaiseError("Failed to write converted coordinate text file.");
    }

    char *line = NULL;
    size_t line_size = 0;
    size_t lineno = 0;
    while(grt_getline(&line, &line_size, in) != -1){
        ++lineno;
        if(grt_is_comment_or_empty_line(line)){
            if(fputs(line, out) == EOF){
                GRTRaiseError("Failed to write converted coordinate text file.");
            }
            continue;
        }

        const char *first_start, *first_end;
        const char *second_start, *second_end;
        real_t first, second;
        parse_text_coordinate_line(
            line, lineno, input_path, &first, &second,
            &first_start, &first_end, &second_start, &second_end);

        real_t out_first, out_second;
        transform_coordinate_pair(
            first, second, Ctrl->C.lon0, Ctrl->C.lat0,
            direction, &out_first, &out_second);
        size_t prefix_len = (size_t)(first_start - line);
        size_t separator_len = (size_t)(second_start - first_end);
        if((fwrite(line, 1, prefix_len, out) != prefix_len) ||
            (fprintf(out, "%.15g", out_first) < 0) ||
            (fwrite(first_end, 1, separator_len, out) != separator_len) ||
            (fprintf(out, "%.15g", out_second) < 0) ||
            (fputs(second_end, out) == EOF)){
            GRTRaiseError("Failed to write converted coordinate text file.");
        }
    }

    if(ferror(in)){
        GRTRaiseError("Failed to read coordinate text file \"%s\".", input_path);
    }
    GRT_SAFE_FREE_PTR(line);
    if(fclose(in) != 0){
        GRTRaiseError("Failed to close coordinate text file \"%s\".", input_path);
    }
    if(fclose(out) != 0){
        GRTRaiseError("Failed to close converted coordinate text file \"%s\".", output_path);
    }
}


/**
 * 执行坐标转换模块
 *
 * @param[in] argc      命令行参数数量
 * @param[in] argv      命令行参数数组
 * @param[in] direction 坐标转换方向
 * @return 模块执行状态
 */
static int coordinate_transform_main(
    int argc, char **argv, GRT_TRANSFORM_DIRECTION direction)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv, direction);

    const char *input_path = Ctrl->G.active ? Ctrl->G.s_path : Ctrl->Q.s_path;
    GRTCheckFileExist(input_path);
    check_output_path(input_path, Ctrl->O.s_path);

    if(Ctrl->Q.active){
        convert_text_file(input_path, Ctrl->O.s_path, Ctrl, direction);
        free_Ctrl(Ctrl);
        return EXIT_SUCCESS;
    }

    // 先按通用二进制格式复制，保留输入 NetCDF 文件的原始格式和内容
    grt_copy_file(input_path, Ctrl->O.s_path);

    // 只在输出副本上读取和修改坐标
    int ncid;
    NC_CHECK(nc_open(Ctrl->O.s_path, NC_WRITE, &ncid));
    GRT_COORD_INFO info;
    inspect_coordinate_layout(ncid, direction, &info);

    convert_output_coordinates(ncid, &info, Ctrl, direction);
    write_reference_and_rename_coordinates(ncid, info.layout, Ctrl, direction);

    NC_CHECK(nc_close(ncid));
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}


/**
 * 执行局部坐标转地理坐标模块
 *
 * @param[in] argc 命令行参数数量
 * @param[in] argv 命令行参数数组
 * @return 模块执行状态
 */
int xy2geo_main(int argc, char **argv)
{
    return coordinate_transform_main(argc, argv, GRT_XY2GEO);
}


/**
 * 执行地理坐标转局部坐标模块
 *
 * @param[in] argc 命令行参数数量
 * @param[in] argv 命令行参数数组
 * @return 模块执行状态
 */
int geo2xy_main(int argc, char **argv)
{
    return coordinate_transform_main(argc, argv, GRT_GEO2XY);
}
