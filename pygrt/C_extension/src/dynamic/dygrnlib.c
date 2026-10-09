/**
 * @file   dygrnlib.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-03
 */

#include <dirent.h>
#include <string.h>
#include <sys/stat.h>

#include "grt/dynamic/dygrnlib.h"
#include "grt/common/checkerror.h"
#include "grt/common/model.h"

/**
 * 识别基本格林函数文件名
 *
 * @param[in]  name  待检查的文件名
 */
static bool gf_filename(const char *name)
{
    // 位移、深度偏导和距离偏导使用相同的基本源型命名
    const char *prefixes[3] = {"", "z", "r"};
    for(int d = 0; d < 3; ++d) {
        GRT_LOOP_ChnlGrid(im, c) {
            if(GRT_SRC_M_ORDERS[im] == 0 && c == 2) {
                continue;
            }
            char expected[32];
            snprintf(expected, sizeof(expected), "%s%s%c.sac", prefixes[d], GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
            if(!strcmp(expected, name)) {
                return true;
            }
        }
    }
    return false;
}

/**
 * 从节点目录按需读取一个基本格林函数的头段
 *
 * @param[in]   path  库节点目录
 * @param[out]  hd    SAC 参考头段
 */
static bool read_header(const char *path, SACHEAD *hd)
{
    DIR *dir = opendir(path);
    if(!dir) {
        return false;
    }

    // 只读一个有效基本格林函数的头段，忽略目录中的附加文件
    struct dirent *entry;
    bool found = false;
    while((entry = readdir(dir))) {
        if(!gf_filename(entry->d_name)) {
            continue;
        }

        char *file = NULL;
        GRT_SAFE_ASPRINTF(&file, "%s/%s", path, entry->d_name);

        SACTRACE *trace = grt_read_SACTRACE(file, true);
        *hd = trace->hd;
        found = true;
        GRT_SAFE_FREE_PTR(file);
        grt_free_SACTRACE(trace);
        break;
    }
    closedir(dir);
    return found;
}

/**
 * 扫描节点几何、起点和震源介质，并检查库内公共时间采样及虚频率
 * @param[in,out]  lib   动态格林函数库
 * @param[in]      path  库节点目录
 * @param[out]     node  节点路径、几何、起点和震源介质
 */
static bool read_node_header(DYGRNLIB *lib, const char *path, DYGRNLIB_NODE *node)
{
    SACHEAD hd;
    if(!read_header(path, &hd)) {
        return false;
    }

    // 参考分量必须有有效的时间采样和起点
    if(hd.npts <= 0 || hd.delta <= 0 || hd.b == SAC_FLOAT_UNDEF) {
        GRTRaiseError("Invalid GF time sampling at %s.", path);
    }

    // 首个节点建立公共采样，后续节点必须使用相同的样本数和间隔
    if(!lib->nt) {
        lib->nt = hd.npts;
        lib->dt = hd.delta;
        lib->wI = GRT_SACHEAD_GET_IMAG_FREQ(&hd);
        if(lib->wI < 0) {
            GRTRaiseError("Invalid GF imaginary frequency at %s.", path);
        }
        lib->stf_decay = exp(-lib->wI * lib->dt);
        lib->refhead = hd;
    } else if(hd.npts != lib->nt || hd.delta != lib->dt) {
        GRTRaiseError("GF nt/dt mismatch at %s.", path);
    }

    // 所有节点必须使用相同虚频率，才能在统一的补偿规则下叠加
    if(GRT_SACHEAD_GET_IMAG_FREQ(&hd) != lib->wI) {
        GRTRaiseError("GF imaginary frequency mismatch at %s.", path);
    }

    // 节点保留几何、起点和震源介质，公共时间采样及虚频率仅在库中保存
    node->path = strdup(path);
    node->zs = hd.evdp;
    node->zr = -hd.stel * 1e-3;
    node->r = hd.dist;
    node->b = hd.b;
    node->vp = GRT_SACHEAD_GET_SRC_VP(&hd);
    node->vs = GRT_SACHEAD_GET_SRC_VS(&hd);
    node->rho = GRT_SACHEAD_GET_SRC_RHO(&hd);
    if(node->vp <= 0 || node->vs < 0 || node->rho <= 0) {
        GRTRaiseError("Invalid source medium in GF header at %s.", path);
    }
    return true;
}

/**
 * 按源深、台深和震中距排序
 *
 * @param[in]  a  第一个库节点
 * @param[in]  b  第二个库节点
 */
static int compare_node(const void *a, const void *b)
{
    const DYGRNLIB_NODE *x = a, *y = b;
    if(x->zs != y->zs) {
        return x->zs < y->zs ? -1 : 1;
    }
    if(x->zr != y->zr) {
        return x->zr < y->zr ? -1 : 1;
    }
    return x->r == y->r ? 0 : x->r < y->r ? -1 : 1;
}

/**
 * 扫描格林函数库，收集节点元数据、公共时间采样、采样轴和原始模型
 *
 * @param[in]  root       格林函数库根目录或单个节点目录
 * @param[in]  calc_upar  是否加载位移空间导数
 */
DYGRNLIB *grt_dygrnlib_load(const char *root, bool calc_upar)
{
    DYGRNLIB *lib = grt_dygrnlib_alloc(calc_upar);

    // 先判断是否直接给定单个节点目录，否则扫描库根目录
    DYGRNLIB_NODE first = {0};
    if(read_node_header(lib, root, &first)) {
        lib->direct_node_input = true;
        lib->nodes = GRT_SAFE_CALLOC(1, sizeof(*lib->nodes));
        lib->nodes[0] = first;
        lib->nnode = 1;
    } else {
        DIR *dir = opendir(root);
        if(!dir) {
            GRTRaiseError("Cannot open Green library %s.", root);
        }

        // 扫描直属子目录，逐个核对目录名与 SAC 头段中的源深、台深和距离
        struct dirent *entry;
        char *modelname = NULL;
        while((entry = readdir(dir))) {
            if(!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, "..")) {
                continue;
            }
            char *path = NULL;
            GRT_SAFE_ASPRINTF(&path, "%s/%s", root, entry->d_name);

            struct stat st;
            // Windows 使用 stat 检查目录，其他平台保持不跟随符号链接
#if _TEST_WHETHER_WIN32_
            if(stat(path, &st) != 0) {
#else
            if(lstat(path, &st) != 0) {
#endif
                GRTRaiseError("Cannot inspect %s.", path);
            }
            if(S_ISDIR(st.st_mode)) {
                // 从目录名末尾提取三个坐标，剩余部分为模型名，可包含下划线
                char *name = strdup(entry->d_name);
                real_t coords[3];
                for(int a = 2; a >= 0; --a) {
                    char *sep = strrchr(name, '_'), *end;
                    if(!sep) {
                        GRTRaiseError("Invalid Green-function node directory name: %s.", path);
                    }
                    coords[a] = strtod(sep + 1, &end);
                    if(end == sep + 1 || *end) {
                        GRTRaiseError("Invalid Green-function node directory name: %s.", path);
                    }
                    *sep = '\0';
                }
                if(!*name) {
                    GRTRaiseError("Invalid Green-function node directory name: %s.", path);
                }

                // 所有子目录必须属于同一个模型，不取决于目录中是否已有 SAC 文件
                if(!modelname) {
                    modelname = strdup(name);
                } else if(strcmp(modelname, name)) {
                    GRTRaiseError("Mixed models in Green-function library %s: %s and %s.", root, modelname, name);
                }
                GRT_SAFE_FREE_PTR(name);

                // 有格林函数的目录才成为库节点，并用容差处理 SAC 单精度坐标舍入
                DYGRNLIB_NODE node = {0};
                if(read_node_header(lib, path, &node)) {
                    if(!GRT_ISCLOSE(node.zs, coords[0]) || !GRT_ISCLOSE(node.zr, coords[1]) || !GRT_ISCLOSE(node.r, coords[2])) {
                        GRTRaiseError("Green-function geometry does not match directory name: %s.", path);
                    }
                    lib->nodes = GRT_SAFE_REALLOC(lib->nodes, (lib->nnode + 1) * sizeof(*lib->nodes));
                    lib->nodes[lib->nnode++] = node;
                }
            }
            GRT_SAFE_FREE_PTR(path);
        }
        closedir(dir);

        // 按统一模型名读取同名模型文件，不将其他普通文件作为模型候选
        if(modelname) {
            char *model_path = NULL;
            GRT_SAFE_ASPRINTF(&model_path, "%s/%s", root, modelname);
            struct stat st;
            if(stat(model_path, &st) != 0 || !S_ISREG(st.st_mode)) {
                GRTRaiseError("Missing model file in Green-function library: %s.", model_path);
            }
            lib->modarr = grt_read_modarr_from_file(model_path, &lib->nlayer, true, false);
            GRT_SAFE_FREE_PTR(model_path);
            GRT_SAFE_FREE_PTR(modelname);
        }
    }
    if(!lib->nnode) {
        GRTRaiseError("No Green functions in %s.", root);
    }

    // 按源深、台深和距离排序，拒绝重复几何，建立后续三维索引顺序
    qsort(lib->nodes, lib->nnode, sizeof(*lib->nodes), compare_node);
    for(size_t i = 1; i < lib->nnode; ++i) {
        if(compare_node(&lib->nodes[i - 1], &lib->nodes[i]) == 0) {
            GRTRaiseError("Duplicate GF geometry: %s and %s.", lib->nodes[i - 1].path, lib->nodes[i].path);
        }
    }

    // 从节点头段收集三个升序采样轴，相同坐标仅保留一次
    real_t **axes[3] = {&lib->depsrcs, &lib->deprcvs, &lib->rs};
    size_t *sizes[3] = {&lib->ndepsrc, &lib->ndeprcv, &lib->nr};
    for(int a = 0; a < 3; ++a) {
        real_t *values = GRT_SAFE_CALLOC(lib->nnode, sizeof(*values));
        for(size_t i = 0; i < lib->nnode; ++i) {
            values[i] = a == 0 ? lib->nodes[i].zs : a == 1 ? lib->nodes[i].zr : lib->nodes[i].r;
        }
        qsort(values, lib->nnode, sizeof(*values), grt_compare_real_t);

        size_t n = 0;
        for(size_t i = 0; i < lib->nnode; ++i) {
            if(!n || values[i] != values[n - 1]) {
                values[n++] = values[i];
            }
        }
        *axes[a] = GRT_SAFE_REALLOC(values, n * sizeof(*values));
        *sizes[a] = n;
    }

    // 节点几何无重复，节点数等于三个轴长度的乘积即可保证所有组合完整，逐次除法避免乘积溢出
    size_t nper_src = lib->nnode / lib->ndepsrc;
    if(lib->nnode % lib->ndepsrc || nper_src % lib->ndeprcv || nper_src / lib->ndeprcv != lib->nr) {
        GRTRaiseError("Incomplete dynamic Green-function library %s: %zu nodes for dimensions %zu x %zu x %zu "
                      "(source depth, receiver depth, distance). "
                      "Every source/receiver depth pair must contain the same distances.",
                      root, lib->nnode, lib->ndepsrc, lib->ndeprcv, lib->nr);
    }
    return lib;
}

/**
 * 求格林函数库内最小正采样间隔
 *
 * @param[in]  lib  格林函数库
 */
real_t grt_dygrnlib_default_subfault_size(const DYGRNLIB *lib)
{
    // 三个采样轴独立取相邻间隔，忽略单点采样轴
    real_t value = INFINITY;
    const real_t *axes[3] = {lib->depsrcs, lib->deprcvs, lib->rs};
    size_t sizes[3] = {lib->ndepsrc, lib->ndeprcv, lib->nr};
    for(int a = 0; a < 3; ++a) {
        for(size_t i = 1; i < sizes[a]; ++i) {
            real_t step = axes[a][i] - axes[a][i - 1];
            if(step > 0 && step < value) {
                value = step;
            }
        }
    }
    if(value == INFINITY) {
        GRTRaiseError("Cannot infer subdivision size; set +idL/dW explicitly.");
    }
    return value;
}

/**
 * 查询三维插值角点及权重
 *
 * @param[in]   lib      格林函数库
 * @param[in]   zs       源深度，km
 * @param[in]   zr       接收深度，km
 * @param[in]   r        震中距，km
 * @param[in]   mode     采样轴查询方式
 * @param[out]  indices  非零权重角点在 nodes 数组中的索引，最多八个
 * @param[out]  weights  角点权重
 */
int grt_dygrnlib_find_corners(
    const DYGRNLIB *lib, real_t zs, real_t zr, real_t r, GRT_SAMPLE_MODE mode,
    size_t indices[8], real_t weights[8])
{
    // 先定位三个采样轴，库加载时已保证所有坐标组合完整
    real_t q[3] = {zs, zr, r}, w[3];
    const real_t *axes[3] = {lib->depsrcs, lib->deprcvs, lib->rs};
    size_t lengths[3] = {lib->ndepsrc, lib->ndeprcv, lib->nr}, idx[3][2];
    for(int a = 0; a < 3; ++a) {
        if(!grt_locate_samples(axes[a], lengths[a], q[a], mode, &idx[a][0], &idx[a][1], &w[a])) {
            return -1;
        }
    }

    // 三个采样轴的左右节点组合成最多八个角点，仅保留非零权重
    int n = 0;
    for(int a = 0; a < 2; ++a) {
        for(int b = 0; b < 2; ++b) {
            for(int c = 0; c < 2; ++c) {
                real_t weight = (a ? w[0] : 1 - w[0]) * (b ? w[1] : 1 - w[1]) * (c ? w[2] : 1 - w[2]);
                if(weight <= 1e-15) {
                    continue;
                }
                size_t index = (idx[0][a] * lib->ndeprcv + idx[1][b]) * lib->nr + idx[2][c];
                indices[n] = index;
                weights[n++] = weight;
            }
        }
    }
    return n;
}

DYGRNLIB_WAVEFORMS *grt_dygrnlib_read_node(const DYGRNLIB *lib, size_t index, const FINITE_FAULT *fault)
{
    const DYGRNLIB_NODE *node = &lib->nodes[index];
    DYGRNLIB_WAVEFORMS *waveforms = GRT_SAFE_CALLOC(1, sizeof(*waveforms));

    // 按位移、深度偏导和距离偏导读取需要的分量，缺失的基本源型无需加载
    const char *prefixes[3] = {"", "z", "r"};
    for(int d = 0; d < (lib->calc_upar ? 3 : 1); ++d) {
        GRT_LOOP_ChnlGrid(im, c) {
            if(!fault->required_gf_sources[im] || (GRT_SRC_M_ORDERS[im] == 0 && c == 2)) {
                continue;
            }
            char *path = NULL;
            GRT_SAFE_ASPRINTF(&path, "%s/%s%s%c.sac", node->path, prefixes[d], GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);

            // 各线程独立拥有当前角点波形，读取和释放均不访问共享状态
            SACTRACE *trace = grt_read_SACTRACE(path, false);

            // 必须核对样本数和间隔，避免按库级 nt 访问长度不足的分量
            if(trace->hd.npts != lib->nt || trace->hd.delta != lib->dt) {
                GRTRaiseError("GF nt/dt mismatch at %s.", path);
            }
            waveforms->gf[d][im][c] = trace->data;
            trace->data = NULL;
            grt_free_SACTRACE(trace);
            GRT_SAFE_FREE_PTR(path);
        }
    }
    return waveforms;
}

void grt_dygrnlib_free_waveforms(DYGRNLIB_WAVEFORMS *waveforms)
{
    if(!waveforms) {
        return;
    }

    // 释放当前角点中实际读取的位移与空间导数，未加载分量保持空指针
    for(int d = 0; d < 3; ++d) {
        GRT_LOOP_ChnlGrid(im, c) {
            GRT_SAFE_FREE_PTR(waveforms->gf[d][im][c]);
        }
    }
    GRT_SAFE_FREE_PTR(waveforms);
}

/**
 * 释放库元数据
 *
 * @param[in,out]  lib  格林函数库
 */
void grt_dygrnlib_free(DYGRNLIB *lib)
{
    if(!lib) {
        return;
    }

    // 库只拥有节点、坐标轴和模型，临时波形由调用方管理
    for(size_t i = 0; i < lib->nnode; ++i) {
        GRT_SAFE_FREE_PTR(lib->nodes[i].path);
    }
    GRT_SAFE_FREE_PTR(lib->nodes);
    GRT_SAFE_FREE_PTR(lib->depsrcs);
    GRT_SAFE_FREE_PTR(lib->deprcvs);
    GRT_SAFE_FREE_PTR(lib->rs);
    GRT_SAFE_FREE_PTR(lib->modarr);
    GRT_SAFE_FREE_PTR(lib);
}

DYGRNLIB *grt_dygrnlib_alloc(bool calc_upar)
{
    DYGRNLIB *lib = GRT_SAFE_CALLOC(1, sizeof(*lib));
    lib->calc_upar = calc_upar;
    return lib;
}
