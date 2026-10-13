/**
 * @file   grt_faultconvert.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-12
 *
 *    将 USGS/WISP 的 PARAM、FSP 动态有限断层转换为 sCoulomb 格式
 *
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    struct {
        bool active;
        char format;            ///< 输入格式，p 为 PARAM，f 为 FSP
        char *path;             ///< 输入断层文件路径
    } F;
    struct {
        bool active;
        real_t lat;             ///< 参考纬度 (degree)
        real_t lon;             ///< 参考经度 (degree)
    } C;
    struct {
        bool active;
        char *path;             ///< 输出 sCoulomb 文件路径
    } O;
    struct {
        bool active;
        size_t index;           ///< 按文件顺序从 1 编号的断层段序号
    } N;
} GRT_MODULE_CTRL;

/** 子断层输入参数，PARAM 记录首个内部积分点，FSP 记录中心 */
typedef struct {
    real_t geographic[3];       ///< 原始 latitude/longitude/depth (degree/degree/km)
    real_t slip;                ///< 总滑动量 (m)
    real_t rake;                ///< 滑动角 (degree)
    real_t strike;              ///< 走向 (degree)
    real_t dip;                 ///< 倾角 (degree)
    real_t delay;               ///< 破裂延迟 (s)
    real_t rise;                ///< 上升时长 (s)
    real_t fall;                ///< 下降时长 (s)
    size_t line;                ///< 输入行号
    real_t rectangle[6];        ///< 输出顶边 east0/north0/east1/north1 与 top/bot
} FAULT_PATCH;

/** 一个断层段及其格式解析状态 */
typedef struct {
    size_t line;                ///< 段头行号
    size_t nx;                  ///< 沿走向子断层数，从段元数据确定
    size_t ny;                  ///< 沿倾向子断层数，从段元数据确定
    size_t expected;            ///< 声明的子断层总数
    size_t count;               ///< 已读取子断层数
    size_t capacity;            ///< FSP 子断层数组的容量，不依赖 Nsbfs 声明
    real_t length;              ///< 子断层沿走向尺寸 (km)
    real_t width;               ///< 子断层沿倾向尺寸 (km)
    real_t plane_length;        ///< FSP 整段沿走向长度 (km)
    real_t plane_width;         ///< FSP 整段沿倾向宽度 (km)
    real_t strike;              ///< FSP 段走向 (degree)
    real_t dip;                 ///< FSP 段倾角 (degree)
    real_t boundary[5][3];       ///< PARAM 边界积分点 latitude/longitude/depth (degree/degree/km)
    size_t nboundary;           ///< 已读取边界积分点数
    bool boundary_header;       ///< 是否读取 PARAM 边界列标题
    bool data_header;           ///< 是否读取子断层列标题
    size_t ncolumns;            ///< FSP 数据列数
    int columns[9];             ///< FSP LAT/LON/X/Y/Z/SLIP/RAKE/TRUP/RISE 所在列
    FAULT_PATCH *patches;        ///< 子断层数组
} FAULT_SEGMENT;


/** 释放参数控制结构体 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    // 释放输入、输出路径，再释放参数结构体
    GRT_SAFE_FREE_PTR(Ctrl->F.path);
    GRT_SAFE_FREE_PTR(Ctrl->O.path);

    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印使用说明 */
static void print_help()
{
// 输出模块语法和选项说明，PARAM 与 FSP 的区别集中在 -F 下
printf("\n[grt faultconvert] %s\n\n", GRT_VERSION);printf(
"    Convert USGS/WISP Kinematic Fault Files to sCoulomb format.\n"
"\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt faultconvert -F[p|f]/<fault> -C<lat>/<lon> -O<file> [-N<index>] [-h]\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -F[p|f]/<fault>\n"
"                  Input USGS/WISP kinematic fault file. Select a format:\n"
"                  + PARAM\n"
"                    set -Fp/<fault>. Recover the first patch center using\n"
"                    half the difference between the third boundary point\n"
"                    and the last data point. Convert slip from cm to m\n"
"                    and preserve rise/fall times.\n"
"                  + FSP\n"
"                    set -Ff/<fault>. Require center coordinates and\n"
"                    asymmetric cosine SVF. Support global and per-segment\n"
"                    sizes. Locate each segment from its first patch center.\n"
"                    Count data rows against the rectangular grid; warn\n"
"                    if only the Nsbfs declaration is inconsistent.\n"
"                    RISE only contains the total duration: use equal\n"
"                    rise/fall times and warn about this approximation.\n"
"                    Other SVFs, top-center and time-window models fail.\n"
"                  Build one undeformed rectangular plane per segment\n"
"                  by projecting its first patch center once, then using\n"
"                  strike/dip and nominal sizes.\n"
"                  Patch order is along strike first, then down dip.\n"
"                  Require uniform strike [0,360] and dip (0,90].\n"
"\n"
"    -C<lat>/<lon>\n"
"                  Reference latitude/longitude (degree), as in geo2xy.\n"
"                  Latitude in (-90,90), longitude in [-180,180].\n"
"\n"
"    -O<file>      Output sCoulomb file, different from the input.\n"
"                  Kode=100, rake/net-slip (m), distances in km,\n"
"                  with a complete -Dc/T1/T2+dTRUP on every row.\n"
"\n"
"    -N<index>     Select one segment by file order, starting at 1.\n"
"                  Omit to convert all segments. Preserve rupture times.\n"
"\n"
"    -h           Display this help message.\n"
"\n"
);
}


/** 从命令行中读取选项并完成基本检查 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    int opt;
    char extra;

    // 逐项解析输入格式、参考坐标、输出路径和断层段选择
    while((opt = getopt(argc, argv, ":F:C:O:N:h")) != -1) {
        switch(opt) {

            // -F[p|f]/<fault>：输入文件及格式
            case 'F':
                if((optarg[0] != 'p' && optarg[0] != 'f') || optarg[1] != '/' || !optarg[2]) {
                    GRTBadOptionError(F, "Expect p/<fault> or f/<fault>.");
                }
                Ctrl->F.active = true;
                Ctrl->F.format = optarg[0];
                GRT_SAFE_FREE_PTR(Ctrl->F.path);
                Ctrl->F.path = strdup(optarg + 2);
                break;

            // -C<lat>/<lon>：局部坐标参考点
            case 'C':
                if(sscanf(optarg, "%lf/%lf %c", &Ctrl->C.lat, &Ctrl->C.lon, &extra) != 2 ||
                   Ctrl->C.lat <= -90.0 || Ctrl->C.lat >= 90.0 || Ctrl->C.lon < -180.0 || Ctrl->C.lon > 180.0) {
                    GRTBadOptionError(C, "Expect latitude/longitude in (-90,90)/[-180,180].");
                }
                Ctrl->C.active = true;
                break;

            // -O<file>：输出 sCoulomb 文件
            case 'O':
                Ctrl->O.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->O.path);
                Ctrl->O.path = strdup(optarg);
                break;

            // -N<index>：仅转换指定断层段，从 1 编号
            case 'N':
                if(!isdigit((unsigned char)optarg[0]) || sscanf(optarg, "%zu%c", &Ctrl->N.index, &extra) != 1 || Ctrl->N.index == 0) {
                    GRTBadOptionError(N, "Segment index must be a positive integer.");
                }
                Ctrl->N.active = true;
                break;

            GRT_Common_Options_in_Switch(optopt);
        }
    }

    // 检查必选项和输入文件，防止输出覆盖原始断层文件
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, F);
    GRTCheckOptionActive(Ctrl, C);
    GRTCheckOptionActive(Ctrl, O);
    GRTCheckFileExist(Ctrl->F.path);
    if(grt_is_same_file(Ctrl->F.path, Ctrl->O.path)) {
        GRTBadOptionError(O, "Output must differ from the input fault file.");
    }
}


/**
 * 从元数据行提取 key = value，避免将 inDx 或 HypX 等字段误认为尺寸
 * @param[in]  Ctrl     参数控制结构体
 * @param[in]  segment  断层段序号
 * @param[in]  line     输入行号
 * @param[in]  text     元数据行
 * @param[in]  key      字段名
 * @param[out] value    数值
 * @return     是否找到字段
 */
static bool read_field(const GRT_MODULE_CTRL *Ctrl, size_t segment, size_t line, const char *text, const char *key, real_t *value)
{
    // 逐次匹配字段名，先排除较长字段名中包含的同名子串
    size_t length = strlen(key);
    for(const char *p = text; (p = strstr(p, key)); p += length) {
        if(p > text && (isalnum((unsigned char)p[-1]) || p[-1] == '_')) {
            continue;
        }

        // 字段名后允许空白，但必须紧接等号
        const char *q = p + length;
        while(isspace((unsigned char)*q)) {
            ++q;
        }
        if(*q++ != '=') {
            continue;
        }

        // 读取等号后的数值，拒绝数值后混入非空白字符
        char *end;
        *value = strtod(q, &end);
        if(end == q || (*end && !isspace((unsigned char)*end))) {
            GRTRaiseError("%s: segment %zu, line %zu: Invalid numeric metadata field.", Ctrl->F.path, segment, line);
        }
        return true;
    }
    return false;
}


/**
 * 检查并读取元数据中的正整数
 * @param[in] Ctrl     参数控制结构体
 * @param[in] segment  断层段序号
 * @param[in] line     输入行号
 * @param[in] value    已读取的数值
 * @return    正整数
 */
static size_t positive_count(const GRT_MODULE_CTRL *Ctrl, size_t segment, size_t line, real_t value)
{
    // 断层段数和剖分数量均要求为正整数
    if(value <= 0.0 || floor(value) != value) {
        GRTRaiseError("%s: segment %zu, line %zu: Expected a positive integer count.", Ctrl->F.path, segment, line);
    }

    return (size_t)value;
}


/**
 * 读取恰好指定数量的数值列，不接受缺列或额外数据
 * @param[in]  Ctrl     参数控制结构体
 * @param[in]  segment  断层段序号
 * @param[in]  line     输入行号
 * @param[in]  text     数据行
 * @param[in]  count    数值列数
 * @param[out] values   数值数组
 */
static void read_values(const GRT_MODULE_CTRL *Ctrl, size_t segment, size_t line, const char *text, size_t count, real_t *values)
{
    // 按指定列数逐个解析，缺列或数值格式错误立即报错
    const char *p = text;
    for(size_t i = 0; i < count; ++i) {
        while(isspace((unsigned char)*p)) {
            ++p;
        }
        char *end;
        values[i] = strtod(p, &end);
        if(end == p || (*end && !isspace((unsigned char)*end))) {
            GRTRaiseError("%s: segment %zu, line %zu: Invalid numeric column or missing rupture parameters.", Ctrl->F.path, segment, line);
        }
        p = end;
    }

    // 读取完指定列后只允许空白，不接受额外的数据列
    while(isspace((unsigned char)*p)) {
        ++p;
    }
    if(*p) {
        GRTRaiseError("%s: segment %zu, line %zu: Unexpected extra data columns.", Ctrl->F.path, segment, line);
    }
}


/**
 * 将经纬度和深度转为内部位置，明确区分投影接口 north/east 和输出 east/north 的顺序
 * @param[in]  Ctrl      参数控制结构体
 * @param[in]  segment   断层段序号
 * @param[in]  line      输入行号
 * @param[in]  lat       纬度 (degree)
 * @param[in]  lon       经度 (degree)
 * @param[in]  depth     深度 (km)
 * @param[out] position  east/north/depth (km)
 */
static void project_position(const GRT_MODULE_CTRL *Ctrl, size_t segment, size_t line,
    real_t lat, real_t lon, real_t depth, real_t position[3])
{
    // 在投影前校验纬度，报错时保留断层段和原始行号
    if(lat < -90.0 || lat > 90.0) {
        GRTRaiseError("%s: segment %zu, line %zu: Latitude must be in [-90,90].", Ctrl->F.path, segment, line);
    }

    // WISP 的几何生成不规范化经度，跨日期变更线时可能超过 180 度，交给公共投影处理
    grt_geo2xy(lat, lon, Ctrl->C.lat, Ctrl->C.lon, &position[1], &position[0]);
    position[2] = depth;
}


/**
 * 在读取下一段或文件结束时检查上一段的结构和尺寸
 * @param[in]     Ctrl     参数控制结构体
 * @param[in]     index    断层段序号
 * @param[in,out] segment  断层段
 */
static void check_segment(const GRT_MODULE_CTRL *Ctrl, size_t index, FAULT_SEGMENT *segment)
{
    // 子断层标题和数量声明必须存在，标称尺寸必须为正
    if(!segment->data_header || segment->expected == 0) {
        GRTRaiseError("%s: segment %zu, line %zu: Missing column header or subfault count declaration.", Ctrl->F.path, index, segment->line);
    }
    if(segment->length <= 0.0 || segment->width <= 0.0) {
        GRTRaiseError("%s: segment %zu, line %zu: Along-strike and down-dip patch sizes must be positive.", Ctrl->F.path, index, segment->line);
    }

    // PARAM 按声明的剖分数量读取，五个边界积分点仍须符合文件结构
    if(Ctrl->F.format == 'p') {
        if(segment->nboundary != 5 || segment->count != segment->expected) {
            GRTRaiseError("%s: segment %zu, line %zu: PARAM requires five boundary points and a complete rectangular grid.",
                Ctrl->F.path, index, segment->line);
        }
    }

    // FSP 同时检查整段尺寸与子断层剖分的一致性
    if(Ctrl->F.format == 'f') {
        if(segment->plane_length <= 0.0 || segment->plane_width <= 0.0) {
            GRTRaiseError("%s: segment %zu, line %zu: FSP segment LEN and WID must be positive.", Ctrl->F.path, index, segment->line);
        }

        // 单段可直接使用 Nx/Nz，多段必须按各段 LEN/WID 与 Dx/Dz 求出整数剖分数
        if(!segment->nx) {
            segment->nx = positive_count(Ctrl, index, segment->line, round(segment->plane_length / segment->length));
        }
        if(!segment->ny) {
            segment->ny = positive_count(Ctrl, index, segment->line, round(segment->plane_width / segment->width));
        }

        // 实际数据行数必须覆盖整个规则网格，不能由错误的 Nsbfs 决定网格数量
        real_t length = segment->nx * segment->length;
        real_t width = segment->ny * segment->width;
        if(fabs(segment->plane_length - length) > 1e-4 * GRT_MAX(1.0, length) ||
           fabs(segment->plane_width - width) > 1e-4 * GRT_MAX(1.0, width) ||
           segment->nx * segment->ny != segment->count) {
            GRTRaiseError("%s: segment %zu, line %zu: FSP sizes and counts must describe one complete rectangular grid.",
                Ctrl->F.path, index, segment->line);
        }

        // 实际行数与网格吻合时，允许历史文件的 Nsbfs 声明错误
        if(segment->expected != segment->count) {
            GRTRaiseWarning("%s: segment %zu, line %zu: Nsbfs declares %zu subfaults, but the complete grid contains %zu; use the actual count.",
                Ctrl->F.path, index, segment->line, segment->expected, segment->count);
        }
    }
}


/**
 * 解析 PARAM 边界或数据列标题，列顺序严格采用 WISP 生成格式
 * @param[in] Ctrl      参数控制结构体
 * @param[in] index     断层段序号
 * @param[in] line      输入行号
 * @param[in] text      列标题
 * @param[in] boundary  是否为边界列标题
 */
static void check_param_columns(const GRT_MODULE_CTRL *Ctrl, size_t index, size_t line, const char *text, bool boundary)
{
    // 边界和子断层使用两套固定列标题
    const char *data_fields[] = {"#Lat.", "Lon.", "depth", "slip", "rake", "strike", "dip", "t_rup", "t_ris", "t_fal", "mo"};
    const char *boundary_fields[] = {"#Lon.", "Lat.", "Depth"};
    const char **expected = boundary ? boundary_fields : data_fields;
    size_t expected_count = boundary ? 3 : 11;

    // 先检查列数，再逐列匹配顺序和名称
    size_t count;
    char **columns = grt_string_split(text, " \t", &count);
    if(count != expected_count) {
        GRTRaiseError("%s: segment %zu, line %zu: Expected the WISP PARAM column header.", Ctrl->F.path, index, line);
    }

    for(size_t i = 0; i < count; ++i) {
        if(strcmp(columns[i], expected[i])) {
            GRTRaiseError("%s: segment %zu, line %zu: Unsupported PARAM column layout.", Ctrl->F.path, index, line);
        }
        free(columns[i]);
    }
    free(columns);
}


/**
 * 读取 PARAM 的分段、边界和子断层记录
 * @param[in]  Ctrl    参数控制结构体
 * @param[in]  lines   输入文本行
 * @param[in]  nlines  输入行数
 * @param[out] count   断层段数
 * @return     新分配的断层段数组
 */
static FAULT_SEGMENT *read_param(const GRT_MODULE_CTRL *Ctrl, char **lines, size_t nlines, size_t *count)
{
    // 保存各段数据，segment 指向当前正在解析的断层段
    FAULT_SEGMENT *segments = NULL;
    FAULT_SEGMENT *segment = NULL;
    size_t total = 0;

    for(size_t i = 0; i < nlines; ++i) {
        // 去掉行首空白并跳过空行，输入行号仍使用原始序号
        const char *text = lines[i];
        while(isspace((unsigned char)*text)) {
            ++text;
        }
        if(!*text) {
            continue;
        }

        real_t value;
        if(strncmp(text, "#Total number of fault_segments", 30) == 0) {
            // 文件头声明总段数，同一文件只接受一次声明
            if(total || !read_field(Ctrl, *count, i + 1, text, "fault_segments", &value)) {
                GRTRaiseError("%s: segment 0, line %zu: Invalid PARAM segment count.", Ctrl->F.path, i + 1);
            }
            total = positive_count(Ctrl, 0, i + 1, value);
        } else if(strncmp(text, "#Fault_segment", 14) == 0) {
            // 开始新段前校验上一段，再分配当前段的存储
            if(segment) {
                check_segment(Ctrl, *count, segment);
            }
            segments = GRT_SAFE_REALLOC(segments, (*count + 1) * sizeof(*segments));
            segment = &segments[(*count)++];
            memset(segment, 0, sizeof(*segment));
            segment->line = i + 1;

            // 从段头读取两个方向的剖分数和子断层尺寸
            size_t id;
            real_t nx, ny;
            int end = 0;
            if(sscanf(text, "#Fault_segment = %zu nx(Along-strike)= %lf Dx = %lfkm ny(downdip)= %lf Dy = %lfkm %n",
                &id, &nx, &segment->length, &ny, &segment->width, &end) != 5 || !end || text[end]) {
                GRTRaiseError("%s: segment %zu, line %zu: Invalid PARAM segment header.", Ctrl->F.path, *count, i + 1);
            }

            segment->nx = positive_count(Ctrl, *count, i + 1, nx);
            segment->ny = positive_count(Ctrl, *count, i + 1, ny);
            segment->expected = segment->nx * segment->ny;
            segment->patches = GRT_SAFE_CALLOC(segment->expected, sizeof(*segment->patches));
        } else if(strncmp(text, "#Lon.", 5) == 0) {
            // 边界列标题必须先于子断层数据列标题
            if(!segment || segment->boundary_header || segment->data_header) {
                GRTRaiseError("%s: segment %zu, line %zu: Unexpected PARAM boundary header.", Ctrl->F.path, *count, i + 1);
            }
            check_param_columns(Ctrl, *count, i + 1, text, true);
            segment->boundary_header = true;
        } else if(strncmp(text, "#Lat.", 5) == 0) {
            // 五个边界点读取完毕后，才允许切换到子断层数据
            if(!segment || segment->nboundary != 5 || segment->data_header) {
                GRTRaiseError("%s: segment %zu, line %zu: Missing boundary points or repeated PARAM column header.", Ctrl->F.path, *count, i + 1);
            }
            check_param_columns(Ctrl, *count, i + 1, text, false);
            segment->data_header = true;
        } else if(*text != '#') {
            // 非注释行按当前状态解释为边界点或子断层记录
            if(!segment || !segment->boundary_header) {
                GRTRaiseError("%s: segment %zu, line %zu: Expected a PARAM segment and boundary header.", Ctrl->F.path, *count, i + 1);
            }

            if(!segment->data_header) {
                if(segment->nboundary == 5) {
                    GRTRaiseError("%s: segment %zu, line %zu: Missing PARAM column header.", Ctrl->F.path, *count, i + 1);
                }

                // PARAM 边界按 longitude/latitude/depth 排列，读入后统一为纬度在前
                real_t values[3];
                read_values(Ctrl, *count, i + 1, text, 3, values);
                real_t *boundary = segment->boundary[segment->nboundary++];
                boundary[0] = values[1];
                boundary[1] = values[0];
                boundary[2] = values[2];
            } else {
                if(segment->count >= segment->expected) {
                    GRTRaiseError("%s: segment %zu, line %zu: PARAM subfault count exceeds the segment declaration.", Ctrl->F.path, *count, i + 1);
                }

                // 逐行保留机制和破裂过程，滑动量由 cm 换算为 m
                real_t values[11];
                read_values(Ctrl, *count, i + 1, text, 11, values);
                FAULT_PATCH *patch = &segment->patches[segment->count++];
                patch->line = i + 1;
                memcpy(patch->geographic, values, sizeof(patch->geographic));
                patch->slip = values[3] * 0.01;
                patch->rake = values[4];
                patch->strike = values[5];
                patch->dip = values[6];
                patch->delay = values[7];
                patch->rise = values[8];
                patch->fall = values[9];
            }
        }
    }

    // 文件结束时检查段数，并校验最后一段的完整性
    if(!segment || total != *count) {
        GRTRaiseError("%s: segment 0, line %zu: Missing PARAM header or inconsistent segment count.", Ctrl->F.path, nlines);
    }
    check_segment(Ctrl, *count, segment);
    return segments;
}


/**
 * 解析 FSP 列标题，并拒绝逐时间窗的动态模型
 * @param[in]     Ctrl     参数控制结构体
 * @param[in]     index    断层段序号
 * @param[in]     line     输入行号
 * @param[in]     text     列标题，不含 % 标记
 * @param[in,out] segment  断层段
 */
static void read_fsp_columns(const GRT_MODULE_CTRL *Ctrl, size_t index, size_t line, const char *text, FAULT_SEGMENT *segment)
{
    // 记录必需字段的位置，允许其它数值列存在以及部分字段调整顺序
    const char *required[] = {"LAT", "LON", "X==EW", "Y==NS", "Z", "SLIP", "RAKE", "TRUP", "RISE"};
    char **columns = grt_string_split(text, " \t", &segment->ncolumns);
    for(size_t k = 0; k < sizeof(required) / sizeof(*required); ++k) {
        segment->columns[k] = -1;
    }

    // 拒绝逐时间窗模型，并防止同一个必需字段重复出现
    for(size_t j = 0; j < segment->ncolumns; ++j) {
        if(strncasecmp(columns[j], "TW", 2) == 0 || strncasecmp(columns[j], "rkTW", 4) == 0) {
            GRTRaiseError("%s: segment %zu, line %zu: Time-window FSP models are not supported.", Ctrl->F.path, index, line);
        }
        for(size_t k = 0; k < sizeof(required) / sizeof(*required); ++k) {
            if(strcmp(columns[j], required[k]) == 0) {
                if(segment->columns[k] >= 0) {
                    GRTRaiseError("%s: segment %zu, line %zu: Repeated FSP data column.", Ctrl->F.path, index, line);
                }
                segment->columns[k] = (int)j;
            }
        }
        free(columns[j]);
    }
    free(columns);

    // 必需字段全部找到后，按实际读取过程扩容，不由 Nsbfs 决定容量
    for(size_t k = 0; k < sizeof(required) / sizeof(*required); ++k) {
        if(segment->columns[k] < 0) {
            GRTRaiseError("%s: segment %zu, line %zu: Missing WISP FSP coordinate, slip or rupture columns.", Ctrl->F.path, index, line);
        }
    }
    segment->capacity = 64;
    segment->patches = GRT_SAFE_CALLOC(segment->capacity, sizeof(*segment->patches));
    segment->data_header = true;
}


/**
 * 读取 FSP，兼容 WISP 历史全局尺寸及当前逐段尺寸布局
 * @param[in]  Ctrl    参数控制结构体
 * @param[in]  lines   输入文本行
 * @param[in]  nlines  输入行数
 * @param[out] count   断层段数
 * @return     新分配的断层段数组
 */
static FAULT_SEGMENT *read_fsp(const GRT_MODULE_CTRL *Ctrl, char **lines, size_t nlines, size_t *count)
{
    // 全局元数据作为单段模型的默认值，多段模型再由各段声明覆盖
    FAULT_SEGMENT *segments = NULL;
    FAULT_SEGMENT *segment = NULL;
    FAULT_SEGMENT global = {.strike = -1.0, .dip = -1.0};
    size_t total = 0;
    bool cosine = false;
    bool centers = false;
    bool source_model = false;

    for(size_t i = 0; i < nlines; ++i) {
        // 去掉行首空白并跳过空行，区分 % 元数据行与数值行
        const char *text = lines[i];
        while(isspace((unsigned char)*text)) {
            ++text;
        }
        if(!*text) {
            continue;
        }

        if(*text == '%') {
            ++text;
            while(isspace((unsigned char)*text)) {
                ++text;
            }

            // 确认源时间函数和坐标语义，不能只根据列名推测输入含义
            real_t value;
            if(strncmp(text, "SVF", 3) == 0) {
                const char *p = strchr(text, ':');
                char name[64];
                if(!p || sscanf(p + 1, "%63s", name) != 1 ||
                   (strcasecmp(name, "Asymetriccosine") && strcasecmp(name, "AsymmetricCosine"))) {
                    GRTRaiseError("%s: segment 0, line %zu: Only the WISP asymmetric cosine SVF is supported.", Ctrl->F.path, i + 1);
                }
                cosine = true;
            }

            if(strncmp(text, "Coordinates are given for", 25) == 0) {
                if(strstr(text, "top-center") || !strstr(text, "center of each")) {
                    GRTRaiseError("%s: segment 0, line %zu: Only WISP center coordinates are supported.", Ctrl->F.path, i + 1);
                }
                centers = true;
            }

            // 源模型标记之后，尺寸和数据标题属于各个断层段
            if(strncmp(text, "SOURCE MODEL PARAMETERS", 23) == 0) {
                source_model = true;
            }

            if(!source_model) {
                // 读取文件级段数、网格尺寸和第一段机制
                if(read_field(Ctrl, 0, i + 1, text, "Nsg", &value)) {
                    total = positive_count(Ctrl, 0, i + 1, value);
                }
                if(read_field(Ctrl, 0, i + 1, text, "Dx", &value)) {
                    global.length = value;
                }
                if(read_field(Ctrl, 0, i + 1, text, "Dz", &value)) {
                    global.width = value;
                }
                if(read_field(Ctrl, 0, i + 1, text, "Nx", &value)) {
                    global.nx = positive_count(Ctrl, 0, i + 1, value);
                }
                if(read_field(Ctrl, 0, i + 1, text, "Nz", &value)) {
                    global.ny = positive_count(Ctrl, 0, i + 1, value);
                }
                if(strncmp(text, "Size", 4) == 0) {
                    read_field(Ctrl, 0, i + 1, text, "LEN", &global.plane_length);
                    read_field(Ctrl, 0, i + 1, text, "WID", &global.plane_width);
                }
                if(strncmp(text, "Mech", 4) == 0) {
                    read_field(Ctrl, 0, i + 1, text, "STRK", &global.strike);
                    read_field(Ctrl, 0, i + 1, text, "DIP", &global.dip);
                }
            } else if(strncmp(text, "SEGMENT #", 9) == 0) {
                // 结束上一段，新段先继承兼容旧版格式所需的全局子断层尺寸
                if(segment) {
                    check_segment(Ctrl, *count, segment);
                }
                segments = GRT_SAFE_REALLOC(segments, (*count + 1) * sizeof(*segments));
                segment = &segments[(*count)++];
                *segment = global;
                segment->line = i + 1;

                // 全局尺寸数量只描述第一段，多段的每一段须使用自己的完整长度和宽度
                segment->nx = segment->ny = 0;
                segment->plane_length = segment->plane_width = 0.0;

                // 多段机制必须由每段声明，不能沿用全局的第一段机制
                segment->strike = segment->dip = -1.0;
                if(!read_field(Ctrl, *count, i + 1, text, "STRIKE", &segment->strike) ||
                   !read_field(Ctrl, *count, i + 1, text, "DIP", &segment->dip)) {
                    GRTRaiseError("%s: segment %zu, line %zu: Missing FSP segment mechanism.", Ctrl->F.path, *count, i + 1);
                }
            } else if(read_field(Ctrl, *count, i + 1, text, "Nsbfs", &value)) {
                // 单段 FSP 没有 SEGMENT 标题，在 Nsbfs 处建立唯一断层段
                if(!segment) {
                    if(total != 1) {
                        GRTRaiseError("%s: segment 0, line %zu: Missing FSP segment headers.", Ctrl->F.path, i + 1);
                    }
                    segments = GRT_SAFE_CALLOC(1, sizeof(*segments));
                    segment = segments;
                    *segment = global;
                    segment->line = i + 1;
                    *count = 1;
                }

                if(segment->expected) {
                    GRTRaiseError("%s: segment %zu, line %zu: Repeated FSP subfault count.", Ctrl->F.path, *count, i + 1);
                }
                segment->expected = positive_count(Ctrl, *count, i + 1, value);
            } else if(segment) {
                // 读取当前段的尺寸，兼容新版逐段给定 Dx/Dz 的布局
                if(read_field(Ctrl, *count, i + 1, text, "Dx", &value)) {
                    segment->length = value;
                }
                if(read_field(Ctrl, *count, i + 1, text, "Dz", &value)) {
                    segment->width = value;
                }
                if(read_field(Ctrl, *count, i + 1, text, "LEN", &value)) {
                    segment->plane_length = value;
                }
                if(read_field(Ctrl, *count, i + 1, text, "WID", &value)) {
                    segment->plane_width = value;
                }

                // 列标题须位于数量声明之后，只允许读取一次
                char first[16], second[16];
                if(sscanf(text, "%15s %15s", first, second) == 2 && strcmp(first, "LAT") == 0 && strcmp(second, "LON") == 0) {
                    if(segment->data_header || !segment->expected) {
                        GRTRaiseError("%s: segment %zu, line %zu: Repeated FSP column header or missing subfault count.",
                            Ctrl->F.path, *count, i + 1);
                    }
                    read_fsp_columns(Ctrl, *count, i + 1, text, segment);
                }
            }
        } else {
            // 数值行必须属于已经声明数量和列布局的段
            if(!segment || !segment->data_header || !segment->expected) {
                GRTRaiseError("%s: segment %zu, line %zu: Expected WISP FSP metadata and data columns.", Ctrl->F.path, *count, i + 1);
            }
            if(segment->count == segment->capacity) {
                segment->capacity *= 2;
                segment->patches = GRT_SAFE_REALLOC(segment->patches, segment->capacity * sizeof(*segment->patches));
            }

            // 根据标题中记录的位置提取坐标、滑动量和破裂时间
            real_t values[segment->ncolumns];
            read_values(Ctrl, *count, i + 1, text, segment->ncolumns, values);
            FAULT_PATCH *patch = &segment->patches[segment->count++];
            patch->line = i + 1;
            patch->geographic[0] = values[segment->columns[0]];
            patch->geographic[1] = values[segment->columns[1]];
            patch->geographic[2] = values[segment->columns[4]];
            patch->slip = values[segment->columns[5]];
            patch->rake = values[segment->columns[6]];
            patch->strike = segment->strike;
            patch->dip = segment->dip;
            patch->delay = values[segment->columns[7]];

            // FSP 仅保留总时长，不能恢复非对称程度，采用等长上升和下降
            patch->rise = patch->fall = 0.5 * values[segment->columns[8]];
        }
    }

    // 检查文件级语义和段数，并校验最后一段
    if(!segment || total != *count || !cosine || !centers || !source_model) {
        GRTRaiseError("%s: segment 0, line %zu: Missing WISP FSP metadata or inconsistent segment count.", Ctrl->F.path, nlines);
    }
    check_segment(Ctrl, *count, segment);
    return segments;
}


/**
 * 恢复首块子断层中心，并仅将该中心投影一次
 *
 * PARAM 第三边界点与末条记录是最后一块的末个、首个内部积分点，差值的一半给出中心偏移
 * WISP 的经纬度转换为线性关系，可直接在 latitude/longitude/depth 中加上该偏移
 * FSP 首条记录已经是中心，两种格式都从首块中心按固定导出顺序展开
 *
 * @param[in]  Ctrl     参数控制结构体
 * @param[in]  index    断层段序号
 * @param[in]  segment  断层段
 * @param[out] position 首块中心，坐标为 east/north/depth
 */
static void project_first_center(const GRT_MODULE_CTRL *Ctrl, size_t index, const FAULT_SEGMENT *segment, real_t position[3])
{
    const FAULT_PATCH *first = &segment->patches[0];
    real_t center[3];
    memcpy(center, first->geographic, sizeof(center));

    // 在经纬度中恢复 PARAM 首块中心，不需要原始定位点或内部积分点数
    if(Ctrl->F.format == 'p') {
        const real_t *corner = segment->boundary[2];
        const real_t *last = segment->patches[segment->count - 1].geographic;
        center[0] += 0.5 * (corner[0] - last[0]);
        center[1] += 0.5 * remainder(corner[1] - last[1], 360.0);
        center[2] += 0.5 * (corner[2] - last[2]);
    }

    project_position(Ctrl, index, first->line, center[0], center[1], center[2], position);
}


/**
 * 将整个矩形等距剖分，并按 WISP 的导出顺序逐行赋予滑动与破裂过程
 *
 * @param[in]     Ctrl     参数控制结构体
 * @param[in]     index    断层段序号
 * @param[in,out] segment  断层段
 */
static void convert_segment(const GRT_MODULE_CTRL *Ctrl, size_t index, FAULT_SEGMENT *segment)
{
    // 矩形段以第一块的走向和倾角构造统一的几何
    const FAULT_PATCH *first = &segment->patches[0];
    if(first->strike < 0.0 || first->strike > 360.0 || first->dip <= 0.0 || first->dip > 90.0) {
        GRTRaiseError("%s: segment %zu, line %zu: Require strike in [0,360] and dip in (0,90].", Ctrl->F.path, index, first->line);
    }

    // 由首块中心确定顶边起点，以及每块沿走向、倾向的增量
    real_t origin[3];
    project_first_center(Ctrl, index, segment, origin);
    real_t sins = sin(first->strike * DEG1);
    real_t coss = cos(first->strike * DEG1);
    real_t sind = sin(first->dip * DEG1);
    real_t cosd = cos(first->dip * DEG1);
    const real_t along[3] = {segment->length * sins, segment->length * coss, 0.0};
    const real_t down[3] = {segment->width * cosd * coss, -segment->width * cosd * sins, segment->width * sind};
    for(int k = 0; k < 3; ++k) {
        origin[k] -= 0.5 * (along[k] + down[k]);
    }

    for(size_t j = 0; j < segment->count; ++j) {
        // 校验各块机制的一致性，以及滑动量和破裂过程的取值
        FAULT_PATCH *patch = &segment->patches[j];
        if(fabs(patch->strike - first->strike) > 1e-6 || fabs(patch->dip - first->dip) > 1e-6) {
            GRTRaiseError("%s: segment %zu, line %zu: Each rectangular segment must have one uniform strike and dip.",
                Ctrl->F.path, index, patch->line);
        }
        if(patch->slip < 0.0 || patch->delay < 0.0 || patch->rise <= 0.0 || patch->fall <= 0.0 ||
           patch->geographic[0] < -90.0 || patch->geographic[0] > 90.0) {
            GRTRaiseError("%s: segment %zu, line %zu: Require valid latitude, nonnegative slip/delay and positive rise/fall durations.",
                Ctrl->F.path, index, patch->line);
        }

        // 将滑动角规范到 [-180,180]，保持等效的滑动方向
        patch->rake = fmod(patch->rake, 360.0);
        if(patch->rake > 180.0) {
            patch->rake -= 360.0;
        }
        if(patch->rake < -180.0) {
            patch->rake += 360.0;
        }

        // WISP 的数组按 (dip, strike) 展平，走向索引变化最快，无需对原始坐标排序
        size_t iL = j % segment->nx;
        size_t iW = j / segment->nx;
        real_t start[3];
        for(int k = 0; k < 3; ++k) {
            start[k] = origin[k] + iL * along[k] + iW * down[k];
        }

        // sCoulomb 用顶边两个端点及上、下边界深度描述矩形
        real_t *r = patch->rectangle;
        r[0] = start[0];
        r[1] = start[1];
        r[2] = start[0] + along[0];
        r[3] = start[1] + along[1];
        r[4] = start[2];
        r[5] = start[2] + down[2];
        if(0.5 * (r[4] + r[5]) < 0.0) {
            GRTRaiseError("%s: segment %zu, line %zu: Subfault center depth must be nonnegative.", Ctrl->F.path, index, patch->line);
        }

        // 输出精度不能使合法矩形退化，或把正时长舍入为零
        if((round(r[0] * 1e6) == round(r[2] * 1e6) && round(r[1] * 1e6) == round(r[3] * 1e6)) ||
           round(r[5] * 1e6) <= round(r[4] * 1e6) || round(patch->rise * 1e6) <= 0.0 || round(patch->fall * 1e6) <= 0.0) {
            GRTRaiseError("%s: segment %zu, line %zu: Geometry or duration is too small for six-decimal output.", Ctrl->F.path, index, patch->line);
        }
    }
}


/** 模块主函数 */
int faultconvert_main(int argc, char **argv)
{
    // 解析命令行并读取原始文件
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);
    FILE *input = GRTCheckOpenFile(Ctrl->F.path, "r");
    size_t nlines;
    char **lines = grt_string_from_file(input, &nlines);
    if(ferror(input)) {
        GRTRaiseError("Failed to read %s.", Ctrl->F.path);
    }
    fclose(input);

    // 按输入格式建立断层段数组，解析完成后释放原始文本
    size_t count = 0;
    FAULT_SEGMENT *segments = Ctrl->F.format == 'p' ? read_param(Ctrl, lines, nlines, &count) : read_fsp(Ctrl, lines, nlines, &count);
    for(size_t i = 0; i <= nlines; ++i) {
        free(lines[i]);
    }
    free(lines);

    // 确定需要输出的段，并在打开输出文件前完成几何和破裂参数转换
    if(Ctrl->N.active && Ctrl->N.index > count) {
        GRTBadOptionError(N, "Segment index exceeds the file's %zu segments.", count);
    }
    size_t begin = Ctrl->N.active ? Ctrl->N.index - 1 : 0;
    size_t end = Ctrl->N.active ? Ctrl->N.index : count;
    for(size_t i = begin; i < end; ++i) {
        convert_segment(Ctrl, i + 1, &segments[i]);
    }
    if(Ctrl->F.format == 'f') {
        GRTRaiseWarning("FSP only preserves total RISE; approximate each patch with equal rise/fall durations (RISE/2).");
    }

    // 所有选中段校验完毕后才打开输出，格式和列宽与 faultdiv 保持一致
    FILE *output = GRTCheckOpenFile(Ctrl->O.path, "w");

    // 写入标准两行表头
    const char *headers[] = {"#", "X-start", "Y-start", "X-fin", "Y-fin", "Kode", "rake", "net-slip", "dip", "top", "bot"};
    for(int row = 0; row < 2; ++row) {
        for(size_t j = 0; j < sizeof(headers) / sizeof(*headers); ++j) {
            int width = j == 0 ? 4 : (j == 5 ? 5 : 15);
            const char *label = row == 0 ? headers[j] : "---------------";
            fprintf(output, "%s%*.*s", j ? " " : "", width, width, label);
        }
        fprintf(output, "\n");
    }

    // 每个子断层写入一条 Kode=100 记录，末尾附上完整的时间函数选项
    for(size_t i = begin; i < end; ++i) {
        FAULT_SEGMENT *segment = &segments[i];
        for(size_t j = 0; j < segment->count; ++j) {
            FAULT_PATCH *patch = &segment->patches[j];
            real_t *r = patch->rectangle;
            fprintf(output, "%4u %15.6f %15.6f %15.6f %15.6f %5u %15.6f %15.6f %15.6f %15.6f %15.6f  -Dc/%.6f/%.6f+d%.6f\n",
                1u, r[0], r[1], r[2], r[3], 100u, patch->rake, patch->slip, patch->dip, r[4], r[5], patch->rise, patch->fall, patch->delay);
        }
    }

    // 确认输出写入成功，再释放段数组及选项结构体
    bool failed = ferror(output);
    if(fclose(output) != 0 || failed) {
        GRTRaiseError("Failed to write %s.", Ctrl->O.path);
    }
    for(size_t i = 0; i < count; ++i) {
        free(segments[i].patches);
    }
    free(segments);
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
