/**
 * @file   syn_output.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-10-04
 *
 *    动态合成结果的 SAC 头段和文件保存
 */

#include <ctype.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>

#include "grt/dynamic/syn_output.h"
#include "grt/dynamic/signals.h"
#include "grt/common/model.h"
#include "grt/common/checkerror.h"

/**
 * 填写接收坐标、机制和有限接收断层的索引头段
 * @param[in]      output  合成输出设置
 * @param[in]      ipt     接收点索引
 * @param[in,out]  hd      SAC 头段
 * @param[in]      depsrc  点源深度，km，有限源为 SAC_FLOAT_UNDEF
 * @param[in]      nlayer  原始模型层数
 * @param[in]      modarr  原始模型矩阵，无模型时为 NULL
 */
void grt_syn_output_set_receiver_header(const DY_SYN_OUTPUT *output, size_t ipt, SACHEAD *hd, real_t depsrc,
                                       size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL])
{
    // SAC 方位角和距离以全局原点为参考，北向和东向坐标另存于接收字段
    const RCV_POINT *rcvs = output->rcvs;
    hd->evdp = depsrc;
    real_t az = atan2(rcvs[ipt].east, rcvs[ipt].north) / DEG1;
    if(az < 0) {
        az += 360;
    }
    hd->az = az;
    hd->baz = fmod(az + 180, 360);
    hd->dist = hypot(rcvs[ipt].north, rcvs[ipt].east);

    // 按实际深度填写介质，界面取上层；无原始模型时保留调用方提供的头段
    if(modarr) {
        real_t top = 0;
        size_t i = 0;
        for(; i + 1 < nlayer; ++i) {
            if(rcvs[ipt].depth <= top + modarr[i][0]) {
                break;
            }
            top += modarr[i][0];
        }
        GRT_SACHEAD_SET_RCV_VP(hd,     modarr[i][1]);
        GRT_SACHEAD_SET_RCV_VS(hd,     modarr[i][2]);
        GRT_SACHEAD_SET_RCV_RHO(hd,    modarr[i][3]);
        GRT_SACHEAD_SET_RCV_QP_INV(hd, modarr[i][4] > 0 ? 1.0 / modarr[i][4] : 0);
        GRT_SACHEAD_SET_RCV_QS_INV(hd, modarr[i][5] > 0 ? 1.0 / modarr[i][5] : 0);
        if(depsrc != SAC_FLOAT_UNDEF) {
            real_t vp, vs, rho;
            grt_modarr_medium_at_depth(nlayer, modarr, depsrc, &vp, &vs, &rho);
            GRT_SACHEAD_SET_SRC_VP(hd,  vp);
            GRT_SACHEAD_SET_SRC_VS(hd,  vs);
            GRT_SACHEAD_SET_SRC_RHO(hd, rho);
        }
    }

    // 坐标和逐点机制使用固定头段字段，缺省字段保持 SAC 未定义值
    GRT_SACHEAD_SET_RCV_NORTH(hd,             rcvs[ipt].north);
    GRT_SACHEAD_SET_RCV_EAST(hd,              rcvs[ipt].east);
    GRT_SACHEAD_SET_RCV_STRIKE(hd,            SAC_FLOAT_UNDEF);
    GRT_SACHEAD_SET_RCV_DIP(hd,               SAC_FLOAT_UNDEF);
    GRT_SACHEAD_SET_RCV_RAKE(hd,              SAC_FLOAT_UNDEF);
    GRT_SACHEAD_SET_RCV_FAULT_INDEX(hd,       SAC_INT_UNDEF);
    hd->stel = -rcvs[ipt].depth * 1e3;

    // 普通接收点仅在输入文件提供机制时填写走向、倾角和滑动角
    if(rcvs[ipt].has_mechanism) {
        GRT_SACHEAD_SET_RCV_STRIKE(hd, rcvs[ipt].strike);
        GRT_SACHEAD_SET_RCV_DIP(hd,    rcvs[ipt].dip);
        GRT_SACHEAD_SET_RCV_RAKE(hd,   rcvs[ipt].rake);
    }

    // 接收点仅记录所属断层的索引
    const FINITE_FAULT *fault = rcvs[ipt].fault;
    if(fault) {
        size_t ifault = fault - output->rcv_faults;
        GRT_SACHEAD_SET_RCV_STRIKE(hd,            fault->strike);
        GRT_SACHEAD_SET_RCV_DIP(hd,               fault->dip);
        GRT_SACHEAD_SET_RCV_RAKE(hd,              fault->rake == GRT_FINITE_FAULT_UNDEFINED_RAKE ? SAC_INT_UNDEF : fault->rake);
        GRT_SACHEAD_SET_RCV_FAULT_INDEX(hd,       ifault);
    }
}

/**
 * 按接收子目录名中的数字索引排序
 * @param[in]  a  第一个目录路径的地址
 * @param[in]  b  第二个目录路径的地址
 */
static int compare_receiver_directories(const void *a, const void *b)
{
    const char *pa = *(const char *const *)a;
    const char *pb = *(const char *const *)b;
    const char *na = strrchr(pa, '/');
    const char *nb = strrchr(pb, '/');
    unsigned long long ia = strtoull(na ? na + 1 : pa, NULL, 10);
    unsigned long long ib = strtoull(nb ? nb + 1 : pb, NULL, 10);
    return (ia > ib) - (ia < ib);
}

char **grt_syn_output_receiver_directories(const char *path, size_t *count)
{
    DIR *dir = opendir(path);
    if(!dir) {
        GRTRaiseError("Cannot open synthesis directory %s.", path);
    }

    // 有子目录时逐个处理子目录，否则直接处理传入的目录
    char **dirs = NULL;
    *count = 0;
    struct dirent *entry;
    while((entry = readdir(dir))) {
        if(entry->d_name[0] == '.') {
            continue;
        }
        char *child = NULL;
        GRT_SAFE_ASPRINTF(&child, "%s/%s", path, entry->d_name);
        struct stat st;
        if(stat(child, &st) == 0 && S_ISDIR(st.st_mode)) {
            dirs = GRT_SAFE_REALLOC(dirs, (*count + 1) * sizeof(*dirs));
            dirs[(*count)++] = child;
        } else {
            GRT_SAFE_FREE_PTR(child);
        }
    }
    closedir(dir);
    if(*count > 0) {
        // 在读取阶段统一恢复接收顺序，供所有动态后处理模块直接使用
        qsort(dirs, *count, sizeof(*dirs), compare_receiver_directories);
    } else {
        dirs = GRT_SAFE_CALLOC(1, sizeof(*dirs));
        dirs[0] = strdup(path);
        *count = 1;
    }
    return dirs;
}

int grt_syn_output_npts(int grn_nt, real_t dt, size_t nfault, const FINITE_FAULT *faults)
{
    int nt = grn_nt;
    for(size_t i = 0; i < nfault; ++i) {
        long long shift = llround(faults[i].stf_delay / dt);
        nt = GRT_MAX(nt, shift + faults[i].stf_npts);
    }
    return nt;
}

void grt_syn_output_save_signal(const char *root, real_t dt, real_t stf_decay,
                                size_t nfault, const FINITE_FAULT *faults, bool finite)
{
    GRTCheckMakeDir(root);

    // 信号覆盖完整破裂过程，每条断层的共享 STF 只累加一次
    int nt = grt_syn_output_npts(1, dt, nfault, faults);
    SACTRACE *sig = grt_new_SACTRACE(dt, nt, 0);
    real_t growth = 1 / stf_decay;
    for(size_t i = 0; i < nfault; ++i) {
        const FINITE_FAULT *fault = &faults[i];
        real_t scale = 1;
        if(finite) {
            scale = 0;
            for(size_t isub = 0; isub < fault->nW * fault->nL; ++isub) {
                scale += fault->scalar_moments[isub];
            }
        }

        // 点源保存归一化 STF，有限源保存完整矩率，写出前恢复物理时间函数
        int shift = (int)llround(fault->stf_delay / dt);
        real_t damping = 1;
        for(int n = 0; n < fault->stf_npts; ++n) {
            sig->data[shift + n] += fault->stfd[n] * damping * scale;
            damping *= growth;
        }
    }

    char *path = NULL;
    GRT_SAFE_ASPRINTF(&path, "%s/sig.sac", root);
    grt_write_SACTRACE(path, sig);
    GRT_SAFE_FREE_PTR(path);
    grt_free_SACTRACE(sig);
}

void grt_syn_output_save_receiver(const DY_SYN_OUTPUT *output, size_t ir, SACTRACE *trace, const real_t *data,
                                  real_t dt, int int_times, int dif_times)
{
    GRTCheckMakeDir(output->root);
    char *outdir = NULL;

    // 极坐标单点直接写根目录，-Q/-U 按接收索引和坐标建立独立子目录
    if(!output->rcv_subdirs) {
        outdir = strdup(output->root);
    } else {
        // 目录名显示两位小数，避免接近零的坐标显示为 -0.00
        const RCV_POINT *receiver = &output->rcvs[ir];
        real_t north = fabs(receiver->north) < 0.005 ? 0 : receiver->north;
        real_t east  = fabs(receiver->east)  < 0.005 ? 0 : receiver->east;
        real_t depth = fabs(receiver->depth) < 0.005 ? 0 : receiver->depth;
        GRT_SAFE_ASPRINTF(&outdir, "%s/%04zu_%.2f_%.2f_%.2f", output->root, ir, north, east, depth);
        GRTCheckMakeDir(outdir);
    }

    // 位移排在前三个分量，其后按求导方向排列九个空间导数
    const char *channels = output->channels;
    int nc = output->calc_upar ? 12 : 3, nt = trace->hd.npts;
    for(int c = 0; c < nc; ++c) {
        // 所有时间算子在工作精度下执行，保存时才转换为 SAC 精度
        memcpy(trace->data, data + (size_t)c * nt, sizeof(*trace->data)*nt);
        for(int j = 0; j < int_times; ++j) {
            grt_trap_integral(trace->data, nt, dt);
        }
        for(int j = 0; j < dif_times; ++j) {
            grt_differential(trace->data, nt, dt);
        }

        // 各分量采用相同的时间算子，空间导数文件名以方向字符为前缀
        char name[8];
        if(c < 3) {
            snprintf(name, sizeof(name), "%c", channels[c]);
        } else {
            snprintf(name, sizeof(name), "%c%c", tolower(channels[(c - 3) / 3]), channels[(c - 3) % 3]);
        }
        snprintf(trace->hd.kcmpnm, sizeof(trace->hd.kcmpnm), "%s", name);

        char *path = NULL;
        GRT_SAFE_ASPRINTF(&path, "%s/%s.sac", outdir, name);
        grt_write_SACTRACE(path, trace);
        GRT_SAFE_FREE_PTR(path);
    }
    GRT_SAFE_FREE_PTR(outdir);
}
