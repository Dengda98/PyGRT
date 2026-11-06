/**
 * @file   grt_synnoise.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-11-01
 *
 *    按固定种子抽样单力源，并分块原位累加长时间 SAC 记录
 */

#include "grt.h"

typedef struct {
    char *id;                   ///< 震源 ID
    real_t north;               ///< 北向坐标，km
    real_t east;                ///< 东向坐标，km
    real_t depth;               ///< 深度，km
    real_t force[3];            ///< 单位力的北、东、向下分量
    real_t strength;            ///< 单力强度，dyne
} NOISE_SOURCE;

typedef struct {
    char *id;                   ///< 台站 ID
    real_t north;               ///< 北向坐标，km
    real_t east;                ///< 东向坐标，km
    real_t depth;               ///< 深度，km
} NOISE_RECEIVER;

typedef struct {
    struct {
        bool active;
        char *path;             ///< 模型路径
        MODEL1D *model;         ///< 原始分层模型
    } M;
    struct {
        bool active;
        char *path;             ///< 震源文件路径
        size_t npoint;          ///< 源点数
        NOISE_SOURCE *sources;  ///< 震源数组
        real_t force_z[2];      ///< 单位力向下分量范围
        real_t azimuth[2];      ///< 从北顺时针的水平角范围，度
        real_t strength[2];     ///< 单力强度范围，dyne
        real_t interval;        ///< 相邻源点的激发间隔，s
        uint64_t seed;          ///< 随机种子
    } C;
    struct {
        bool active;
        char *path;                ///< 台站文件路径
        size_t npoint;             ///< 台站数
        NOISE_RECEIVER *receivers; ///< 台站数组
        real_t dt;                ///< 采样间隔，s
        real_t length;            ///< 完整记录长度，s
        int nt;                   ///< 最终 SAC 样本数
    } U;
    struct {
        bool active;
        char *path;             ///< 输出目录
    } O;
    struct {
        bool active;
        char *option;           ///< 传给 modsum 的完整阶数选项
        ssize_t max_order;      ///< 最大阶数，-1 表示全部
    } N;
    struct {
        bool active;
        int nprocs;             ///< syn 并发进程数，也是 eigenv 和 modsum 的线程数
    } P;
    struct {
        bool active;
    } s;
} GRT_MODULE_CTRL;

/** 释放输入点集、模型和选项 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl)
{
    for(size_t i = 0; i < Ctrl->C.npoint; ++i) {
        GRT_SAFE_FREE_PTR(Ctrl->C.sources[i].id);
    }
    for(size_t i = 0; i < Ctrl->U.npoint; ++i) {
        GRT_SAFE_FREE_PTR(Ctrl->U.receivers[i].id);
    }
    GRT_SAFE_FREE_PTR(Ctrl->C.sources);
    GRT_SAFE_FREE_PTR(Ctrl->U.receivers);

    grt_free_mod1d(Ctrl->M.model);
    GRT_SAFE_FREE_PTR(Ctrl->M.path);
    GRT_SAFE_FREE_PTR(Ctrl->C.path);
    GRT_SAFE_FREE_PTR(Ctrl->U.path);
    GRT_SAFE_FREE_PTR(Ctrl->O.path);
    GRT_SAFE_FREE_PTR(Ctrl->N.option);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(void)
{
printf("\n"
"[grt synnoise] %s\n\n", GRT_VERSION);printf(
"    Simulate ZNE ambient-noise velocity records (cm/s) by modal summation.\n"
"    Each source fires once in file order, with random force direction and strength.\n"
"\n"
"\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt synnoise -M<model> -C<sfile> -U<rfile> -O<dir>\n"
"                 [-N[<n1>[/<n2>[/<dn>]]]] [-P<nprocs>] [-s] [-h]\n"
"\n"
"\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -M<model>    Filepath to 1D horizontally layered halfspace \n"
"                 model. The model file has 6 columns: \n"
"\n"
"         +-------+----------+----------+-------------+----+----+\n"
"         | H(km) | Vp(km/s) | Vs(km/s) | Rho(g/cm^3) | Qp | Qa |\n"
"         +-------+----------+----------+-------------+----+----+\n"
"\n"
"         The 1st column H is the layer thickness. \n"
"         For compatibility, if the \"thickness\" of the first layer \n"
"         set 0.0, then the first column will be the layer top depth, \n"
"         which must be in an ascending order.\n"
"         The number of layers are unlimited.\n"
"\n"
"    -C<sfile>        Source file: sampling parameters followed by source points.\n"
"                     See the input format below.\n"
"\n"
"    -U<rfile>        Receiver file: recording parameters followed by receiver points.\n"
"                     See the input format below.\n"
"\n"
"    -O<dir>          New or empty output directory. SAC files are saved as\n"
"                     <dir>/noise/<receiver_ID>/{Z,N,E}.sac, together with\n"
"                     the model, phase{R,L}.nc, group{R,L}.nc and source/receiver logs.\n"
"                     SAC headers unused6 and unused7 store receiver North and East coordinates (km).\n"
"                     The elevation header stel stores -1000 * depth (m), with depth in km.\n"
"\n"
"    -N[n1[/n2[/dn]]] Mode selection. Default -N0; bare -N selects all modes.\n"
"\n"
"    -P<nprocs>       Maximum concurrent processes. Default: all available cores.\n"
"                     Used as the thread count for eigenv and modsum; each syn uses one thread.\n"
"\n"
"    -s               Suppress regular progress messages.\n"
"\n"
"    -h               Display this help message.\n"
"\n"
"\n"
"Input files:\n"
"----------------------------------------------------------------\n"
"    Source file (-C): the five parameter lines are required in this order.\n"
"    Ranges use min max; equal limits give a fixed value.\n"
"\n"
"        FORCE_Z        -1 1          # Unit-force downward component, within [-1, 1]\n"
"        FORCE_AZIMUTH  0 360         # Degrees clockwise from North, within [0, 360]\n"
"        STRENGTH       1e15 2e15     # Nonnegative force strength (dyne)\n"
"        INTERVAL       1.08          # Positive firing interval (s), starting at 0\n"
"        SEED           20261031      # Nonnegative integer seed\n"
"\n"
"        # ID      north  east  depth\n"
"        SRC001       10     0      0\n"
"        SRC002        0    12      0\n"
"\n"
"    Receiver file (-U): the two parameter lines are required in this order.\n"
"    DT and LENGTH must be positive; LENGTH/DT must be an integer.\n"
"\n"
"        DT      0.1                 # Sampling interval (s)\n"
"        LENGTH  3600                # Record length (s)\n"
"\n"
"        # ID      north  east  depth\n"
"        STA001        0     0      0\n"
"        STA002        0     1    0.5\n"
"\n"
"    Coordinates are in km; depths are nonnegative. IDs must be unique in each file.\n"
"    Blank lines and # comments are accepted. Events at or after length are skipped.\n"
"    All paths, including the grt executable, allow ASCII letters, digits and\n"
"    _ . - / : only (Windows also accepts backslashes).\n"
"\n"
"\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt synnoise -Mmilrow -Csources.txt -Ureceivers.txt -Onoise\n"
"\n");
}

/** 从命令行读取模块选项 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv)
{
    Ctrl->N.option = strdup("-N0");
    int opt;
    char extra;
    while((opt = getopt(argc, argv, ":M:C:U:O:N::P:sh")) != -1) {
        switch(opt) {
            // -M<model>
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.path = strdup(optarg);
                break;

            // -C<sfile>
            case 'C':
                Ctrl->C.active = true;
                Ctrl->C.path = strdup(optarg);
                break;

            // -U<rfile>
            case 'U':
                Ctrl->U.active = true;
                Ctrl->U.path = strdup(optarg);
                break;

            // -O<dir>
            case 'O':
                Ctrl->O.active = true;
                Ctrl->O.path = strdup(optarg);
                break;

            // -N[n1[/n2[/dn]]]
            case 'N':
                Ctrl->N.active = true;
                GRT_SAFE_FREE_PTR(Ctrl->N.option);
                GRT_SAFE_ASPRINTF(&Ctrl->N.option, "-N%s", optarg ? optarg : "");
                Ctrl->N.max_order = -1;
                // 仅提取 eigenv 所需的最大阶数，完整阶数选项交给 modsum 解析
                if(optarg != NULL) {
                    const char *upper = strchr(optarg, '/');
                    if(!grt_string_composed_of(optarg, GRT_NUM_STR "+-/") ||
                       sscanf(upper ? upper + 1 : optarg, "%zd", &Ctrl->N.max_order) != 1) {
                        GRTBadOptionError(N, "");
                    }
                }
                break;

            // -P<nprocs>
            case 'P':
                Ctrl->P.active = true;
                if(sscanf(optarg, "%d%c", &Ctrl->P.nprocs, &extra) != 1 || Ctrl->P.nprocs <= 0) {
                    GRTBadOptionError(P, "Need a positive process count.");
                }
                break;

            // -s
            case 's':
                Ctrl->s.active = true;
                break;

            GRT_Common_Options_in_Switch(optopt);
        }
    }
    GRTCheckOptionActive(Ctrl, M);
    GRTCheckOptionActive(Ctrl, C);
    GRTCheckOptionActive(Ctrl, U);
    GRTCheckOptionActive(Ctrl, O);

    // 命令直接拼接，不支持路径中的空白、引号或 shell 特殊字符
    const char *paths[] = {GRT_EXECUTABLE, Ctrl->M.path, Ctrl->C.path, Ctrl->U.path, Ctrl->O.path};
    const char *allowed = _TEST_WHETHER_WIN32_ ?
                          "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" GRT_NUM_STR "_.-/:\\" :
                          "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" GRT_NUM_STR "_.-/:";
    for(size_t i = 0; i < sizeof(paths) / sizeof(*paths); ++i) {
        if(!*paths[i] || !grt_string_composed_of(paths[i], allowed)) {
            GRTRaiseError("Unsupported path characters: %s. Use ASCII letters, digits, _, ., -, / and : (Windows also accepts backslashes).", paths[i]);
        }
    }
}

/** 按 ID 排序指针，检查重复时不改变输入事件顺序 */
static int compare_noise_id(const void *a, const void *b)
{
    const char *x = *(const char *const *)a;
    const char *y = *(const char *const *)b;
    return strcmp(x, y);
}

/** 按震源深度排序指针，建库分组时不改变事件及其累加顺序 */
static int compare_noise_depth(const void *a, const void *b)
{
    const NOISE_SOURCE *x = *(const NOISE_SOURCE *const *)a;
    const NOISE_SOURCE *y = *(const NOISE_SOURCE *const *)b;
    return (x->depth > y->depth) - (x->depth < y->depth);
}

/**
 * 读取震源或台站文件的参数头和四列坐标
 * @param[in,out]  Ctrl    模块参数
 * @param[in]      source  是否读取震源而非接收点
 */
static void read_noise_points(GRT_MODULE_CTRL *Ctrl, bool source)
{
    const char *path = source ? Ctrl->C.path : Ctrl->U.path;
    FILE *fp = GRTCheckOpenFile(path, "r");
    char *line = NULL;
    size_t capacity = 0, lineno = 0, count = 0, allocated = 0;
    const char *source_parameters[] = {"FORCE_Z", "FORCE_AZIMUTH", "STRENGTH", "INTERVAL", "SEED"};
    const char *receiver_parameters[] = {"DT", "LENGTH"};
    const char **parameters = source ? source_parameters : receiver_parameters;
    size_t nheader = source ? 5 : 2;
    real_t *ranges[] = {Ctrl->C.force_z, Ctrl->C.azimuth, Ctrl->C.strength};
    size_t nparameter = 0;

    while(grt_getline(&line, &capacity, fp) != -1) {
        ++lineno;

        // 去掉行尾注释，空行与纯注释行不参与参数和点集计数
        char *comment = strchr(line, '#');
        if(comment) {
            *comment = '\0';
        }
        if(grt_is_comment_or_empty_line(line)) {
            continue;
        }

        // 震源文件前五行、台站文件前两行是参数头，其余有效行使用相同的四列格式
        char extra;
        if(nparameter < nheader) {
            char name[32];
            if(sscanf(line, "%31s", name) != 1 || strcmp(name, parameters[nparameter])) {
                GRTRaiseError("%s:%zu: expected %s in the %s header.", path, lineno, parameters[nparameter], source ? "source" : "receiver");
            }

            if(!source) {
                real_t *value = nparameter == 0 ? &Ctrl->U.dt : &Ctrl->U.length;
                if(sscanf(line, "%*s %lf %c", value, &extra) != 1 || *value <= 0) {
                    GRTRaiseError("%s:%zu: %s must be positive seconds.", path, lineno, parameters[nparameter]);
                }
            } else if(nparameter < 3) {
                real_t *range = ranges[nparameter];
                if(sscanf(line, "%*s %lf %lf %c", &range[0], &range[1], &extra) != 2 || range[0] > range[1]) {
                    GRTRaiseError("%s:%zu: invalid %s min max range.", path, lineno, parameters[nparameter]);
                }
            } else if(nparameter == 3) {
                if(sscanf(line, "%*s %lf %c", &Ctrl->C.interval, &extra) != 1 || Ctrl->C.interval <= 0) {
                    GRTRaiseError("%s:%zu: INTERVAL must be positive seconds.", path, lineno);
                }
            } else {
                // 种子允许取到 UINT64_MAX，单独按无符号十进制读取
                char value[32];
                if(sscanf(line, "%*s %31s %c", value, &extra) != 1 || !grt_string_composed_of(value, GRT_NUM_STR)) {
                    GRTRaiseError("%s:%zu: SEED must be an unsigned 64-bit decimal integer.", path, lineno);
                }
                errno = 0;
                Ctrl->C.seed = strtoull(value, NULL, 10);
                if(errno == ERANGE) {
                    GRTRaiseError("%s:%zu: SEED must be an unsigned 64-bit decimal integer.", path, lineno);
                }
            }
            ++nparameter;
            continue;
        }

        // ID 用于输出目录名，限制长度并排除目录分隔符和相对路径
        char id[256];
        NOISE_RECEIVER point = {0};
        if(sscanf(line, "%255s %lf %lf %lf %c", id, &point.north, &point.east, &point.depth, &extra) != 4 || point.depth < 0) {
            GRTRaiseError("%s:%zu: expected ID north east nonnegative_depth.", path, lineno);
        }
        if(strchr(id, '/') || strchr(id, '\\') || !strcmp(id, ".") || !strcmp(id, "..")) {
            GRTRaiseError("%s:%zu: ID cannot contain directory separators or be '.' or '..'.", path, lineno);
        }

        // 点集按需扩容，保留文件顺序作为后续的激发顺序
        if(count == allocated) {
            allocated = allocated ? 2 * allocated : 256;
            if(source) {
                Ctrl->C.sources = GRT_SAFE_REALLOC(Ctrl->C.sources, allocated * sizeof(*Ctrl->C.sources));
            } else {
                Ctrl->U.receivers = GRT_SAFE_REALLOC(Ctrl->U.receivers, allocated * sizeof(*Ctrl->U.receivers));
            }
        }
        point.id = strdup(id);
        if(source) {
            Ctrl->C.sources[count] = (NOISE_SOURCE){.id = point.id, .north = point.north, .east = point.east, .depth = point.depth};
        } else {
            Ctrl->U.receivers[count] = point;
        }
        ++count;
    }

    if(ferror(fp) || count == 0) {
        GRTRaiseError("Incomplete or empty point file: %s.", path);
    }
    fclose(fp);
    GRT_SAFE_FREE_PTR(line);

    // 范围已检查升序，此处只检查各物理量的取值边界
    if(source) {
        if(Ctrl->C.force_z[0] < -1 || Ctrl->C.force_z[1] > 1) {
            GRTRaiseError("%s: FORCE_Z must be within [-1, 1].", path);
        }
        if(Ctrl->C.azimuth[0] < 0 || Ctrl->C.azimuth[1] > 360) {
            GRTRaiseError("%s: FORCE_AZIMUTH must be within [0, 360].", path);
        }
        if(Ctrl->C.strength[0] < 0) {
            GRTRaiseError("%s: STRENGTH must be nonnegative.", path);
        }
    } else {
        real_t count = Ctrl->U.length / Ctrl->U.dt;
        if(fabs(count - round(count)) > 1e-6 || count < 1 || count > 2147483647.0) {
            GRTRaiseError("%s: LENGTH/DT must be an integer supported by SAC v6 npts.", path);
        }
        Ctrl->U.nt = (int)llround(count);
    }

    // 只排序指针来检查重复 ID，不改变点集顺序
    const char **sorted = GRT_SAFE_MALLOC(count * sizeof(*sorted));
    for(size_t i = 0; i < count; ++i) {
        sorted[i] = source ? Ctrl->C.sources[i].id : Ctrl->U.receivers[i].id;
    }
    qsort(sorted, count, sizeof(*sorted), compare_noise_id);
    for(size_t i = 1; i < count; ++i) {
        if(!strcmp(sorted[i - 1], sorted[i])) {
            GRTRaiseError("%s: duplicate ID %s.", path, sorted[i]);
        }
    }
    GRT_SAFE_FREE_PTR(sorted);

    if(source) {
        Ctrl->C.npoint = count;
    } else {
        Ctrl->U.npoint = count;
    }
}

/** 使用局部 SplitMix64 状态返回 [0, 1) 均匀随机数，不依赖系统 rand 或内部模块的随机状态 */
static real_t noise_uniform(uint64_t *state)
{
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return (real_t)((value ^ (value >> 31)) >> 11) * 0x1.0p-53;
}

/**
 * 排序并写出去重后的坐标列表
 * @param[in]      path    输出文件路径
 * @param[in,out]  values  待排序的坐标数组
 * @param[in]      count   坐标数量
 * @return                去重后的坐标数量
 */
static size_t write_unique_values(const char *path, real_t *values, size_t count)
{
    qsort(values, count, sizeof(*values), grt_compare_real_t);
    FILE *fp = GRTCheckOpenFile(path, "w");

    // 相邻坐标相同时只写一次，避免生成重复的格林函数节点
    size_t unique = 0;
    for(size_t i = 0; i < count; ++i) {
        if(i == 0 || values[i] != values[i - 1]) {
            fprintf(fp, "%.17g\n", values[i]);
            ++unique;
        }
    }

    if(fclose(fp) != 0) {
        GRTRaiseError("Cannot write geometry file: %s.", path);
    }
    return unique;
}

/** 只递归清理本次创建的临时目录，不通过 shell 执行路径 */
static void remove_noise_tree(const char *path)
{
    DIR *dir = opendir(path);
    if(!dir) {
        GRTRaiseError("Cannot open temporary directory: %s.", path);
    }

    // 先递归删除子目录和文件，再移除当前空目录
    struct dirent *entry;
    while((entry = readdir(dir))) {
        if(!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, "..")) {
            continue;
        }
        char *child = NULL;
        GRT_SAFE_ASPRINTF(&child, "%s/%s", path, entry->d_name);
        struct stat st;
        if(stat(child, &st) != 0) {
            GRTRaiseError("Cannot inspect temporary path: %s.", child);
        }
        if(S_ISDIR(st.st_mode)) {
            remove_noise_tree(child);
        } else if(remove(child) != 0) {
            GRTRaiseError("Cannot remove temporary file: %s.", child);
        }
        GRT_SAFE_FREE_PTR(child);
    }

    closedir(dir);
    if(rmdir(path) != 0) {
        GRTRaiseError("Cannot remove temporary directory: %s.", path);
    }
}

/**
 * 检查有效事件的介质和源台距，并估计格林函数时窗
 * @param[in]   Ctrl    模块参数
 * @param[out]  nevent  记录时间内的事件数
 * @return              格林函数时窗，s
 */
static real_t get_noise_time_window(const GRT_MODULE_CTRL *Ctrl, size_t *nevent)
{
    const MODEL1D *model = Ctrl->M.model;
    if(model->topbound != GRT_BOUND_FREE || model->botbound != GRT_BOUND_HALFSPACE) {
        GRTRaiseError("synnoise requires a free surface and half-space bottom.");
    }

    // 取模型中的最小正 S 波速度，用于估计最慢传播时间
    real_t vmin = INFINITY;
    for(size_t i = 0; i < model->n; ++i) {
        if(model->Vb[i] > 0) {
            vmin = GRT_MIN(vmin, model->Vb[i]);
        }
    }

    // 只检查记录时间内的源，统计有效事件数和最远水平源台距
    real_t rmax = 0;
    size_t count = 0;
    for(size_t is = 0; is < Ctrl->C.npoint; ++is) {
        real_t vp, vs, rho;
        const NOISE_SOURCE *src = &Ctrl->C.sources[is];
        real_t time = is * Ctrl->C.interval;
        if(time >= Ctrl->U.length) {
            GRTRaiseWarning("Skipping %zu sources starting with %s at %.8g s, outside [0, %.8g).",
                            Ctrl->C.npoint - is, src->id, time, Ctrl->U.length);
            break;
        }

        // 面波单力源要求固体介质，且不能与台站具有相同水平坐标
        grt_modarr_medium_at_depth(model->nmodarr, model->modarr, src->depth, &vp, &vs, &rho);
        if(vs == 0) {
            GRTRaiseError("Source %s is in a liquid layer; modal force sources require a solid layer.", src->id);
        }
        for(size_t ir = 0; ir < Ctrl->U.npoint; ++ir) {
            const NOISE_RECEIVER *rcv = &Ctrl->U.receivers[ir];
            real_t distance = hypot(rcv->north - src->north, rcv->east - src->east);
            if(GRT_IS_ZERO(distance)) {
                GRTRaiseError("Zero horizontal distance between source %s and station %s.", src->id, rcv->id);
            }
            rmax = GRT_MAX(rmax, distance);
        }
        ++count;
    }

    // 按两倍最远传播时间估计格林函数时窗，补足末端采样点并向上取偶数
    int npts = GRT_MAX(4, (int)ceil(2 * rmax / (vmin * Ctrl->U.dt)) + 1);
    npts += npts % 2;
    real_t length = npts * Ctrl->U.dt;
    if(!Ctrl->s.active) {
        GRTRaiseInfo("GF: Rmax=%.8g km, vmin=%.8g km/s, npts=%d, length=%.8g s, df=%.8g Hz.", rmax, vmin, npts, length, 1 / length);
    }

    *nevent = count;
    return length;
}

/**
 * 创建输出目录、输入副本，并报告完整 SAC 记录的预计体积
 * @param[in]  Ctrl  模块参数
 */
static void initialize_noise_output(const GRT_MODULE_CTRL *Ctrl)
{
    const MODEL1D *model = Ctrl->M.model;

    // 输出目录必须为空，避免覆盖已有记录
    GRTCheckMakeDir(Ctrl->O.path);
    DIR *out = opendir(Ctrl->O.path);
    if(!out) {
        GRTRaiseError("Cannot open output directory: %s.", Ctrl->O.path);
    }
    struct dirent *entry;
    while((entry = readdir(out))) {
        if(strcmp(entry->d_name, ".") && strcmp(entry->d_name, "..")) {
            GRTRaiseError("Output directory must be empty: %s.", Ctrl->O.path);
        }
    }
    closedir(out);

    // 模型副本不能与固定输出文件重名
    const char *reserved[] = {"noise", ".synnoise_tmp", "sources.info", "receivers.info",
                              "phaseR.nc", "phaseL.nc", "groupR.nc", "groupL.nc"};
    for(size_t i = 0; i < sizeof(reserved) / sizeof(*reserved); ++i) {
        if(!strcmp(model->modelname, reserved[i])) {
            GRTRaiseError("Model filename conflicts with synnoise output: %s.", model->modelname);
        }
    }

    // 复制模型和台站文件，并建立波形输出根目录
    char *model_copy = NULL, *receiver_copy = NULL, *noise_root = NULL;
    GRT_SAFE_ASPRINTF(&model_copy, "%s/%s", Ctrl->O.path, model->modelname);
    GRT_SAFE_ASPRINTF(&receiver_copy, "%s/receivers.info", Ctrl->O.path);
    GRT_SAFE_ASPRINTF(&noise_root, "%s/noise", Ctrl->O.path);
    grt_copy_file(Ctrl->M.path, model_copy);
    grt_copy_file(Ctrl->U.path, receiver_copy);
    GRTCheckMakeDir(noise_root);
    GRT_SAFE_FREE_PTR(model_copy);
    GRT_SAFE_FREE_PTR(receiver_copy);

    // 在逐台站创建长记录前报告预计文件大小
    unsigned long long file_bytes = GRT_SAC_HEADER_BYTES + (unsigned long long)SAC_DATA_SIZEOF * Ctrl->U.nt;
    unsigned long long nfile = 3 * Ctrl->U.npoint;
    unsigned long long total_bytes = nfile * file_bytes;
    if(!Ctrl->s.active) {
        GRTRaiseInfo("Noise SAC: %d samples, %.8g s, %llu files, %llu bytes/file, %llu bytes total (%.3f MiB, %.3f GiB).",
                     Ctrl->U.nt, Ctrl->U.length, nfile, file_bytes, total_bytes, total_bytes / 1048576.0, total_bytes / 1073741824.0);
    }

    GRT_SAFE_FREE_PTR(noise_root);
}

/**
 * 为当前台站分块创建三分量空 SAC 记录
 * @param[in]  Ctrl  模块参数
 * @param[in]  rcv   当前台站
 */
static void initialize_noise_station(const GRT_MODULE_CTRL *Ctrl, const NOISE_RECEIVER *rcv)
{
    char *station_dir = NULL;
    GRT_SAFE_ASPRINTF(&station_dir, "%s/noise/%s", Ctrl->O.path, rcv->id);
    GRTCheckMakeDir(station_dir);

    // 各台站共用模拟零时刻，位置和分量写入 SAC 头段
    SACHEAD hd = new_sac_head(Ctrl->U.dt, Ctrl->U.nt, 0);
    hd.o = 0;
    hd.iztype = IO;
    hd.stel = -1000 * rcv->depth;
    snprintf(hd.kstnm, sizeof(hd.kstnm), "%.8s", rcv->id);
    GRT_SACHEAD_SET_RCV_NORTH(&hd, rcv->north);
    GRT_SACHEAD_SET_RCV_EAST(&hd, rcv->east);

    for(int c = 0; c < 3; ++c) {
        char *path = NULL;
        snprintf(hd.kcmpnm, sizeof(hd.kcmpnm), "%c", "ZNE"[c]);
        hd.cmpinc = c == 0 ? 0 : 90;
        hd.cmpaz = c == 2 ? 90 : 0;
        GRT_SAFE_ASPRINTF(&path, "%s/%c.sac", station_dir, "ZNE"[c]);
        grt_create_zero_SAC(path, &hd);
        GRT_SAFE_FREE_PTR(path);
    }
    GRT_SAFE_FREE_PTR(station_dir);
}

/**
 * 计算 Rayleigh 和 Love 波频散，供后续各源点共用
 * @param[in]  Ctrl    模块参数
 * @param[in]  length  格林函数时窗，s
 */
static void compute_noise_dispersion(const GRT_MODULE_CTRL *Ctrl, real_t length)
{
    char *modes = NULL;
    if(Ctrl->N.max_order < 0) {
        modes = strdup("-N");
    } else {
        GRT_SAFE_ASPRINTF(&modes, "-N%zd", Ctrl->N.max_order);
    }

    // 两种波型共用均匀频率网格，频散计算使用指定的线程数
    int nthreads = Ctrl->P.active ? Ctrl->P.nprocs : grt_get_num_threads(-1);
    for(int w = 0; w < 2; ++w) {
        char *command = NULL;
        GRT_SAFE_ASPRINTF(&command, "%s eigenv -M%s -S%c -F0/%.17g/%.17g -C%s/phase%c.nc %s -P%d -s",
                          GRT_EXECUTABLE, Ctrl->M.path, "RL"[w], 0.5 / Ctrl->U.dt, 1 / length,
                          Ctrl->O.path, "RL"[w], modes, nthreads);
        grt_run_command(command);
        GRT_SAFE_FREE_PTR(command);

        // 计算相速度文件中全部频率、模态对应的群速度频散
        GRT_SAFE_ASPRINTF(&command, "%s eigenfn -C%s/phase%c.nc -N -U%s/group%c.nc",
                          GRT_EXECUTABLE, Ctrl->O.path, "RL"[w], Ctrl->O.path, "RL"[w]);
        grt_run_command(command);
        GRT_SAFE_FREE_PTR(command);
    }
    GRT_SAFE_FREE_PTR(modes);
}

/**
 * 按固定种子抽样全部有效事件，并写出对齐的震源日志
 * @param[in,out]  Ctrl    模块参数及待抽样震源
 * @param[in]      nevent  有效事件数
 */
static void sample_noise_sources(GRT_MODULE_CTRL *Ctrl, size_t nevent)
{
    char *events_path = NULL;
    GRT_SAFE_ASPRINTF(&events_path, "%s/sources.info", Ctrl->O.path);
    FILE *events = GRTCheckOpenFile(events_path, "w");

    // 预先统计坐标、时刻和 ID 的显示宽度，保证表头及所有数据行逐列对齐
    const char *event_columns[] = {"ID", "north", "east", "depth", "time", "force_N/force_E/force_Z", "strength"};
    int event_widths[7];
    for(int i = 0; i < 7; ++i) {
        event_widths[i] = strlen(event_columns[i]);
    }
    event_widths[5] = GRT_MAX(event_widths[5], snprintf(NULL, 0, "%+.6f/%+.6f/%+.6f", 1.0, 1.0, 1.0));
    for(size_t is = 0; is < nevent; ++is) {
        const NOISE_SOURCE *src = &Ctrl->C.sources[is];
        real_t time = is * Ctrl->C.interval;
        int widths[] = {strlen(src->id), snprintf(NULL, 0, "%.6f", src->north), snprintf(NULL, 0, "%.6f", src->east),
                        snprintf(NULL, 0, "%.6f", src->depth), snprintf(NULL, 0, "%.6f", time)};
        for(int i = 0; i < 5; ++i) {
            event_widths[i] = GRT_MAX(event_widths[i], widths[i]);
        }
    }
    for(int i = 0; i < 2; ++i) {
        event_widths[6] = GRT_MAX(event_widths[6], snprintf(NULL, 0, "%.6e", Ctrl->C.strength[i]));
    }
    fprintf(events, "%-*s    %*s    %*s    %*s     %*s     %*s     %*s\n",
            event_widths[0], event_columns[0], event_widths[1], event_columns[1], event_widths[2], event_columns[2],
            event_widths[3], event_columns[3], event_widths[4], event_columns[4], event_widths[5], event_columns[5],
            event_widths[6], event_columns[6]);

    // 按文件顺序串行抽样，每个事件固定消耗三个随机数，使结果不受并发数影响
    uint64_t random_state = Ctrl->C.seed;
    for(size_t is = 0; is < nevent; ++is) {
        NOISE_SOURCE *src = &Ctrl->C.sources[is];
        real_t time = is * Ctrl->C.interval;
        real_t force_z = Ctrl->C.force_z[0] + (Ctrl->C.force_z[1] - Ctrl->C.force_z[0]) * noise_uniform(&random_state);
        real_t angle = (Ctrl->C.azimuth[0] + (Ctrl->C.azimuth[1] - Ctrl->C.azimuth[0]) * noise_uniform(&random_state)) * DEG1;
        real_t strength = Ctrl->C.strength[0] + (Ctrl->C.strength[1] - Ctrl->C.strength[0]) * noise_uniform(&random_state);
        real_t horizontal = sqrt(GRT_MAX(0, 1 - force_z * force_z));

        src->force[0] = horizontal * cos(angle);
        src->force[1] = horizontal * sin(angle);
        src->force[2] = force_z;
        src->strength = strength;

        // 日志仅作显示，合成仍使用未舍入的方向和强度
        char force[64];
        snprintf(force, sizeof(force), "%+.6f/%+.6f/%+.6f", src->force[0], src->force[1], src->force[2]);
        fprintf(events, "%-*s    %*.6f    %*.6f    %*.6f     %*.6f     %*s     %*.6e\n",
                event_widths[0], src->id, event_widths[1], src->north, event_widths[2], src->east,
                event_widths[3], src->depth, event_widths[4], time, event_widths[5], force, event_widths[6], strength);
    }

    if(fclose(events) != 0) {
        GRTRaiseError("Cannot close event table: %s.", events_path);
    }

    GRT_SAFE_FREE_PTR(events_path);
}

/**
 * 将一个源点的三分量短响应平移并累加到台站长记录
 * @param[in]  Ctrl         模块参数
 * @param[in]  receiver_id  台站 ID
 * @param[in]  syn_path     短响应目录
 * @param[in]  time         发震时刻，s
 */
static void add_noise_response(const GRT_MODULE_CTRL *Ctrl, const char *receiver_id, const char *syn_path, real_t time)
{
    long long sample = llround(time / Ctrl->U.dt);

    for(int c = 0; c < 3; ++c) {
        char *response_path = NULL, *target_path = NULL;
        GRT_SAFE_ASPRINTF(&response_path, "%s/%c.sac", syn_path, "ZNE"[c]);
        GRT_SAFE_ASPRINTF(&target_path, "%s/noise/%s/%c.sac", Ctrl->O.path, receiver_id, "ZNE"[c]);

        // 短响应可整体读入，长记录只更新与本次响应重叠的片段
        SACHEAD hd;
        float *response = read_sac(response_path, &hd);
        if(!response || fabs(hd.delta - Ctrl->U.dt) > 1e-6 * Ctrl->U.dt) {
            GRTRaiseError("Invalid sampling in short response: %s.", response_path);
        }

        // 将发震时刻对齐到采样点，并裁剪记录范围以外的响应
        long long first = sample + llround(hd.b / Ctrl->U.dt);
        long long begin = GRT_MAX(0, first);
        long long end = GRT_MIN(Ctrl->U.nt, first + hd.npts);
        if(end > begin) {
            grt_add_SAC_window(target_path, begin, end - begin, response + (begin - first));
        }

        GRT_SAFE_FREE_PTR(response);
        GRT_SAFE_FREE_PTR(response_path);
        GRT_SAFE_FREE_PTR(target_path);
    }
}

/**
 * 逐台站批量建库，再并发合成各震源的短响应并按源文件顺序累加
 * @param[in]  Ctrl     模块参数
 * @param[in]  nevent   有效事件数
 */
static void compute_noise_records(const GRT_MODULE_CTRL *Ctrl, size_t nevent)
{
    char *temporary = NULL, *green = NULL, *dist_path = NULL;
    GRT_SAFE_ASPRINTF(&temporary, "%s/.synnoise_tmp", Ctrl->O.path);
    GRT_SAFE_ASPRINTF(&green, "%s/green", temporary);
    GRT_SAFE_ASPRINTF(&dist_path, "%s/distances.txt", temporary);
    GRTCheckMakeDir(temporary);

    // 仅按原始震源深度分组，避免将深度与距离的笛卡尔积全部计算
    NOISE_SOURCE **sources = GRT_SAFE_MALLOC(nevent * sizeof(*sources));
    real_t *distances = GRT_SAFE_MALLOC(nevent * sizeof(*distances));
    for(size_t is = 0; is < nevent; ++is) {
        sources[is] = &Ctrl->C.sources[is];
    }
    qsort(sources, nevent, sizeof(*sources), compare_noise_depth);

    int nthreads = Ctrl->P.active ? Ctrl->P.nprocs : grt_get_num_threads(-1);
    size_t nprocs = GRT_MIN((size_t)nthreads, nevent);
    char **syn_paths = GRT_SAFE_CALLOC(nprocs, sizeof(*syn_paths));
    GRT_PROCESS *processes = GRT_SAFE_MALLOC(nprocs * sizeof(*processes));
    for(size_t j = 0; j < nprocs; ++j) {
        GRT_SAFE_ASPRINTF(&syn_paths[j], "%s/syn%zu", temporary, j);
    }
    if(!Ctrl->s.active) {
        GRTRaiseInfo("Processing %zu events at %zu stations: modsum uses %d thread(s), syn uses up to %zu concurrent process(es).",
                     nevent, Ctrl->U.npoint, nthreads, nprocs);
    }

    // 台站作为距离计算的中心，源深、台站深和合成时的源台方向均保持原义
    for(size_t ir = 0; ir < Ctrl->U.npoint; ++ir) {
        const NOISE_RECEIVER *rcv = &Ctrl->U.receivers[ir];
        initialize_noise_station(Ctrl, rcv);
        if(!Ctrl->s.active) {
            GRTRaiseInfo("Station %s (%zu/%zu): initialized three empty SAC files.", rcv->id, ir + 1, Ctrl->U.npoint);
        }

        // 每个深度组一次传入当前台站到全部组内震源的距离，并去重重复节点
        for(size_t start = 0; start < nevent;) {
            size_t end = start + 1;
            real_t depth = sources[start]->depth;
            while(end < nevent && sources[end]->depth == depth) {
                ++end;
            }
            for(size_t is = start; is < end; ++is) {
                distances[is - start] = hypot(rcv->north - sources[is]->north, rcv->east - sources[is]->east);
            }
            size_t ndist = write_unique_values(dist_path, distances, end - start);

            // modsum 内部按频率使用 OpenMP，两个波型依次计算以避免线程超额占用
            for(int w = 0; w < 2; ++w) {
                char *command = NULL;
                GRT_SAFE_ASPRINTF(&command, "%s modsum -C%s/phase%c.nc -D%.17g/%.17g -R%s -O%s %s -Gvh -P%d",
                                  GRT_EXECUTABLE, Ctrl->O.path, "RL"[w], depth, rcv->depth, dist_path, green, Ctrl->N.option, nthreads);
                if(!Ctrl->s.active) {
                    GRTRaiseInfo("Station %s: modsum %c GFs (depsrc=%.17g km, ndist=%zu)...", rcv->id, "RL"[w], depth, ndist);
                }
                grt_run_command(command);
                if(!Ctrl->s.active) {
                    GRTRaiseInfo("Station %s: modsum %c done (depsrc=%.17g km).", rcv->id, "RL"[w], depth);
                }
                GRT_SAFE_FREE_PTR(command);
            }
            start = end;
        }

        char progress_prefix[128];
        if(!Ctrl->s.active) {
            snprintf(progress_prefix, sizeof(progress_prefix), "Synthesizing Noise Records (%zu/%zu): ", ir + 1, Ctrl->U.npoint);
            grt_printprogressBar(progress_prefix, 0);
        }
        for(size_t start = 0; start < nevent; start += nprocs) {
            size_t count = GRT_MIN(nprocs, nevent - start);

            // 同一台站的多个震源并发合成，每个子进程使用独立输出目录和一个线程
            for(size_t j = 0; j < count; ++j) {
                const NOISE_SOURCE *src = &Ctrl->C.sources[start + j];
                real_t north = rcv->north - src->north;
                real_t east = rcv->east - src->east;
                real_t distance = hypot(north, east);
                real_t azimuth = atan2(east, north) / DEG1;
                if(azimuth < 0) {
                    azimuth += 360;
                }

                // 只利用震中距的对称性批量建库，合成时仍使用原始震源、台站和方位角
                char *command = NULL;
                GRT_SAFE_ASPRINTF(&command,
                                  "%s syn -G%s/%s_%.17g_%.17g_%.17g -A%.17g -S%.17g -F%.17g/%.17g/%.17g -O%s -N -P1 -s", // 震源为阶跃型，台站记录速度，-I1 和 -J1 抵消
                                  GRT_EXECUTABLE, green, Ctrl->M.model->modelname, src->depth, rcv->depth, distance,
                                  azimuth, src->strength, src->force[0], src->force[1], src->force[2], syn_paths[j]);
                processes[j] = grt_start_command(command);
                GRT_SAFE_FREE_PTR(command);
            }
            grt_wait_processes(count, processes);

            // 主进程按文件顺序更新当前台站的长记录，保持逐次浮点累加顺序不变
            for(size_t j = 0; j < count; ++j) {
                real_t time = (start + j) * Ctrl->C.interval;
                add_noise_response(Ctrl, rcv->id, syn_paths[j], time);
                remove_noise_tree(syn_paths[j]);
            }
            if(!Ctrl->s.active) {
                grt_printprogressBar(progress_prefix, (start + count) * 100 / nevent);
            }
        }

        // 当前台站完成后立即释放格林函数，临时目录和进程槽位供下一台站复用
        remove_noise_tree(green);
        if(!Ctrl->s.active) {
            GRTRaiseInfo("Station %s (%zu/%zu) completed.", rcv->id, ir + 1, Ctrl->U.npoint);
        }
    }

    remove_noise_tree(temporary);
    for(size_t j = 0; j < nprocs; ++j) {
        GRT_SAFE_FREE_PTR(syn_paths[j]);
    }
    GRT_SAFE_FREE_PTR(syn_paths);
    GRT_SAFE_FREE_PTR(processes);
    GRT_SAFE_FREE_PTR(sources);
    GRT_SAFE_FREE_PTR(distances);
    GRT_SAFE_FREE_PTR(dist_path);
    GRT_SAFE_FREE_PTR(green);
    GRT_SAFE_FREE_PTR(temporary);
}

/** 子模块主函数 */
int synnoise_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = GRT_SAFE_CALLOC(1, sizeof(*Ctrl));
    getopt_from_command(Ctrl, argc, argv);

    // 读入震源、台站和模型，在创建输出前完成介质及几何检查
    read_noise_points(Ctrl, true);
    read_noise_points(Ctrl, false);
    Ctrl->M.model = grt_read_mod1d_from_file(Ctrl->M.path, -1, -1, false, false);

    size_t nevent = 0;
    real_t length = get_noise_time_window(Ctrl, &nevent);

    // 初始化输出目录，再计算所有源点共用的频散结果
    initialize_noise_output(Ctrl);
    compute_noise_dispersion(Ctrl, length);

    // 串行抽样后并发合成，保持随机数消耗和记录累加顺序固定
    sample_noise_sources(Ctrl, nevent);
    compute_noise_records(Ctrl, nevent);

    if(!Ctrl->s.active) {
        GRTRaiseInfo("Completed %zu events at %zu stations.", nevent, Ctrl->U.npoint);
    }
    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
