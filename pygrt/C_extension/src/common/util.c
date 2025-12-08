/**
 * @file   util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 其它辅助函数
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "grt/common/util.h"
#include "grt/common/model.h"
#include "grt/common/const.h"
#include "grt/common/sacio.h"
#include "grt/common/myfftw.h"
#include "grt/common/travt.h"

#include "grt/common/checkerror.h"

char ** grt_string_split(const char *string, const char *delim, size_t *size)
{
    char *str_copy = strdup(string);  // 创建字符串副本，以免修改原始字符串
    char *token = strtok(str_copy, delim);

    char **s_split = NULL;
    *size = 0;

    while(token != NULL){
        s_split = (char**)realloc(s_split, sizeof(char*)*(*size+1));
        s_split[*size] = NULL;
        s_split[*size] = (char*)realloc(s_split[*size], sizeof(char)*(strlen(token)+1));
        strcpy(s_split[*size], token);

        token = strtok(NULL, delim);
        (*size)++;
    }
    free(str_copy);

    return s_split;
}

char ** grt_string_from_file(FILE *fp, size_t *size){
    char **s_split = NULL;
    *size = 0;
    s_split = (char**)realloc(s_split, sizeof(char*)*(*size+1));
    s_split[*size] = NULL;

    size_t len=0;
    while(grt_getline(&s_split[*size], &len, fp) != -1){
        s_split[*size][strlen(s_split[*size])-1] = '\0';  // 换行符换为终止符
        (*size)++;
        s_split = (char**)realloc(s_split, sizeof(char*)*(*size+1));
        s_split[*size] = NULL;
    }
    return s_split;
}

bool grt_string_composed_of(const char *str, const char *alws){
    bool allowed[256] = {false};  // 初始全为false（不允许）

    // 标记允许的字符
    for (int i = 0; alws[i] != '\0'; i++) {
        unsigned char c = alws[i];  // 转为无符号避免负数索引
        allowed[c] = true;
    }

    // 检查目标字符串中的每个字符
    for (int i = 0; str[i] != '\0'; i++) {
        unsigned char c = str[i];
        if (!allowed[c]) {  // 若字符不在允许集合中
            return false;
        }
    }

    // 所有字符均在允许集合中
    return true;
}

int grt_string_ncols(const char *string, const char* delim){
    int count = 0;
    
    const char *str = string;
    while (*str) {
        // 跳过所有分隔符
        while (*str && strchr(delim, *str)) str++;
        // 如果还有非分隔符字符，增加计数
        if (*str) count++;
        // 跳过所有非分隔符字符
        while (*str && !strchr(delim, *str)) str++;
    }
    
    return count;
}


const char* grt_get_basename(const char* path) {
    // 找到最后一个 '/'
    char* last_slash = strrchr(path, '/'); 
    
#ifdef _WIN32
    char* last_backslash = strrchr(path, '\\');
    if (last_backslash && (!last_slash || last_backslash > last_slash)) {
        last_slash = last_backslash;
    }
#endif
    if (last_slash) {
        // 返回最后一个 '/' 之后的部分
        return last_slash + 1; 
    }
    // 如果没有 '/'，整个路径就是最后一项
    return path; 
}


void grt_trim_whitespace(char* str) {
    char* end;
    
    // 去除首部空白
    while (isspace((unsigned char)*str)) str++;
    
    // 去除尾部空白
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    
    // 写入终止符
    *(end + 1) = '\0';
}


bool grt_is_comment_or_empty(const char* line) {
    // 跳过前导空白
    while (isspace((unsigned char)*line)) line++;
    
    // 检查是否为空行或注释行
    return (*line == '\0' || *line == GRT_COMMENT_HEAD);
}


ssize_t grt_getline(char **lineptr, size_t *n, FILE *stream){
    if (!lineptr || !n || !stream) {
        return -1;
    }
    
    char *buf = *lineptr;
    size_t size = *n;
    size_t len = 0;
    int c;
    
    // 如果缓冲区为空，分配初始缓冲区
    if (buf == NULL || size == 0) {
        size = 128;
        buf = malloc(size);
        if (buf == NULL) {
            return -1;
        }
    }
    
    // 逐字符读取直到换行符或EOF
    while ((c = fgetc(stream)) != EOF) {
        // 检查是否需要扩展缓冲区
        if (len + 1 >= size) {
            size_t new_size = size * 2;
            char *new_buf = realloc(buf, new_size);
            if (new_buf == NULL) {
                free(buf);
                return -1;
            }
            buf = new_buf;
            size = new_size;
        }
        
        buf[len++] = c;
        
        // 遇到换行符停止读取
        if (c == '\n') {
            break;
        }
    }
    
    // 如果没有读取到任何字符且遇到EOF
    if (len == 0 && c == EOF) {
        return -1;
    }
    
    // 添加字符串终止符
    buf[len] = '\0';
    
    *lineptr = buf;
    *n = size;
    
    return len;
}




/**
 * 将一条数据反变换回时间域再进行处理，保存到SAC文件
 * 
 * @param[in]         srcname       震源类型
 * @param[in]         ch            三分量类型（Z,R,T）
 * @param[in]         delayT        延迟时间
 * @param[in]         wI            虚频率
 * @param[in,out]     pt_fh         FFTW结构体
 * @param[in,out]     pt_hd         SAC头段变量结构体指针
 * @param[in]         s_output_subdir    保存路径所在文件夹
 * @param[in]         s_prefix           sac文件名以及通道名名称前缀
 * @param[in]         sgn                数据待乘符号(-1/1)
 * @param[in]         grncplx   复数形式的格林函数频谱
 * 
 */
static void write_one_to_sac(
    const char *srcname, const char ch, const real_t delayT,
    const real_t wI, GRT_FFTW_HOLDER *pt_fh,
    SACHEAD *pt_hd, const char *s_output_subdir, const char *s_prefix,
    const int sgn, const cplx_t *grncplx)
{
    snprintf(pt_hd->kcmpnm, sizeof(pt_hd->kcmpnm), "%s%s%c", s_prefix, srcname, ch);
    
    char *s_outpath = NULL;
    GRT_SAFE_ASPRINTF(&s_outpath, "%s/%s.sac", s_output_subdir, pt_hd->kcmpnm);

    // 执行fft任务会修改数组，需重新置零
    grt_reset_fftw_holder_zero(pt_fh);
    
    // 赋值复数，包括时移
    cplx_t cfac, ccoef;
    cfac = exp(I*PI2*pt_fh->df*delayT);
    ccoef = sgn;
    // 只赋值有效长度，其余均为0
    for(size_t i=0; i<pt_fh->nf_valid; ++i){
        pt_fh->W_f[i] = grncplx[i] * ccoef;
        ccoef *= cfac;
    }


    if(! pt_fh->naive_inv){
        // 发起fft任务 
        fftw_execute(pt_fh->plan);
    } else {
        grt_naive_inverse_transform_double(pt_fh);
    }
    

    // 归一化，并处理虚频
    // 并转为 SAC 需要的单精度类型
    float *float_arr = (float*)malloc(sizeof(float)*pt_fh->nt);
    real_t fac, coef;
    coef = pt_fh->df * exp(delayT*wI);
    fac = exp(wI*pt_fh->dt);
    for(size_t i=0; i<pt_fh->nt; ++i){
        float_arr[i] = pt_fh->w_t[i] * coef;
        coef *= fac;
    }

    // 以sac文件保存到本地
    write_sac(s_outpath, *pt_hd, float_arr);

    GRT_SAFE_FREE_PTR(float_arr);
    GRT_SAFE_FREE_PTR(s_outpath);
}


/**
 * 处理单个震中距对应的数据逆变换和SAC保存
 * 
 * @param[in]         command       模块名
 * @param[in]         mod1d         模型结构体指针
 * @param[in]         s_prefix      保存路径前缀
 * @param[in]         wI            虚频率
 * @param[in,out]     pt_fh         FFTW结构体
 * @param[in]         s_dist        输入的震中距字符串
 * @param[in]         dist          震中距
 * @param[in]         depsrc        震源深度
 * @param[in]         deprcv        接收深度
 * @param[in]         delayT0       延迟时间
 * @param[in]         delayV0       参考速度
 * @param[in]         calc_upar     是否计算位移偏导
 * @param[in]         doEX          是否保存爆炸源结果
 * @param[in]         doVF          是否保存垂直力源结果
 * @param[in]         doHF          是否保存水平力源结果
 * @param[in]         doDC          是否保存剪切力源结果
 * @param[in]         doDC          是否保存剪切力源结果
 * @param[in,out]     pt_hd         SAC头段变量结构体指针
 * @param[in]         chalst        要保存的分量字符串
 * @param[in]         grn           格林函数频谱结果
 * @param[in]         grn_uiz       格林函数对z偏导频谱结果
 * @param[in]         grn_uir       格林函数对r偏导频谱结果
 * 
 */
static void single_freq2time_write_to_file(
    const char *command, const GRT_MODEL1D *mod1d, const char *s_prefix, 
    const real_t wI, GRT_FFTW_HOLDER *pt_fh,
    const char *s_dist, const real_t dist,
    const real_t depsrc, const real_t deprcv,
    const real_t delayT0, const real_t delayV0, const bool calc_upar,
    const bool doEX, const bool doVF, const bool doHF, const bool doDC, 
    SACHEAD *pt_hd, const char *chalst,
    pt_cplxChnlGrid grn, 
    pt_cplxChnlGrid grn_uiz, 
    pt_cplxChnlGrid grn_uir)
{
    // 文件保存子目录
    char *s_output_subdir = NULL;
    
    GRT_SAFE_ASPRINTF(&s_output_subdir, "%s_%s", s_prefix, s_dist);
    GRTCheckMakeDir(command, s_output_subdir);

    // 时间延迟 
    real_t delayT = delayT0;
    if(delayV0 > 0.0)   delayT += sqrt( GRT_SQUARE(dist) + GRT_SQUARE(deprcv - depsrc) ) / delayV0;
    // 修改SAC头段时间变量
    pt_hd->b = delayT;

    // 计算理论走时
    pt_hd->t0 = grt_compute_travt1d(mod1d->Thk, mod1d->Va, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
    strcpy(pt_hd->kt0, "P");
    pt_hd->t1 = grt_compute_travt1d(mod1d->Thk, mod1d->Vb, mod1d->n, mod1d->isrc, mod1d->ircv, dist);
    strcpy(pt_hd->kt1, "S");

    GRT_LOOP_ChnlGrid(im, c){
        if(! doEX  && im==0)  continue;
        if(! doVF  && im==1)  continue;
        if(! doHF  && im==2)  continue;
        if(! doDC  && im>=3)  continue;

        int modr = GRT_SRC_M_ORDERS[im];
        int sgn=1;  // 用于反转Z分量

        if(modr==0 && GRT_ZRT_CODES[c]=='T')  continue;  // 跳过输出0阶的T分量

        // 判断是否为所需的分量
        if(strchr(chalst, GRT_ZRT_CODES[c]) == NULL)  continue;

        // Z分量反转
        sgn = (GRT_ZRT_CODES[c]=='Z') ? -1 : 1;

        write_one_to_sac(
            GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c], delayT, 
            wI, pt_fh,
            pt_hd, s_output_subdir, "", sgn, 
            grn[im][c]);

        if(calc_upar){
            write_one_to_sac(
                GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c], delayT, 
                wI, pt_fh,
                pt_hd, s_output_subdir, "z", sgn*(-1), 
                grn_uiz[im][c]);

            write_one_to_sac(
                GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c], delayT, 
                wI, pt_fh,
                pt_hd, s_output_subdir, "r", sgn, 
                grn_uir[im][c]);
        }
    }


    GRT_SAFE_FREE_PTR(s_output_subdir);
}



void grt_GF_freq2time_write_to_file(
    const char *command, const GRT_MODEL1D *mod1d, 
    const char *s_output_dir, const char *s_modelname, const char *s_depsrc, const char *s_deprcv,    
    const real_t wI, GRT_FFTW_HOLDER *pt_fh,
    const size_t nr, char *s_dists[nr], const real_t dists[nr], real_t travtPS[nr][2],
    const real_t depsrc, const real_t deprcv,
    const real_t delayT0, const real_t delayV0, const bool calc_upar,
    const bool doEX, const bool doVF, const bool doHF, const bool doDC, 
    const char *chalst,
    pt_cplxChnlGrid grn[nr], 
    pt_cplxChnlGrid grn_uiz[nr], 
    pt_cplxChnlGrid grn_uir[nr])
{
    // 建立SAC头文件，包含必要的头变量
    SACHEAD hd = new_sac_head(pt_fh->dt, pt_fh->nt, delayT0);
    // 发震时刻作为参考时刻
    hd.o = 0.0; 
    hd.iztype = IO; 
    // 记录震源和台站深度
    hd.evdp = depsrc; // km
    hd.stel = (-1.0)*deprcv*1e3; // m
    // 写入虚频率
    hd.user0 = wI;
    // 写入接受点的Vp,Vs,rho
    hd.user1 = mod1d->Va[mod1d->ircv];
    hd.user2 = mod1d->Vb[mod1d->ircv];
    hd.user3 = mod1d->Rho[mod1d->ircv];
    hd.user4 = mod1d->Qainv[mod1d->ircv];
    hd.user5 = mod1d->Qbinv[mod1d->ircv];
    // 写入震源点的Vp,Vs,rho
    hd.user6 = mod1d->Va[mod1d->isrc];
    hd.user7 = mod1d->Vb[mod1d->isrc];
    hd.user8 = mod1d->Rho[mod1d->isrc];

    char *s_output_dirprefx = NULL;
    GRT_SAFE_ASPRINTF(&s_output_dirprefx, "%s/%s_%s_%s", s_output_dir, s_modelname, s_depsrc, s_deprcv);
    
    // 做反傅里叶变换，保存SAC文件
    for(size_t ir=0; ir<nr; ++ir){
        hd.dist = dists[ir];

        single_freq2time_write_to_file(
            command, mod1d, s_output_dirprefx, 
            wI, pt_fh,
            s_dists[ir], dists[ir], depsrc, deprcv,
            delayT0, delayV0, calc_upar,
            doEX, doVF, doHF, doDC,
            &hd, chalst, grn[ir], grn_uiz[ir], grn_uir[ir]);
        
        // 记录走时
        if(travtPS != NULL){
            travtPS[ir][0] = hd.t0;
            travtPS[ir][1] = hd.t1;
        }
    }

    // 输出警告：当震源位于液体层中时，仅允许计算爆炸源对应的格林函数
    if(mod1d->Vb[mod1d->isrc]==0.0){
        GRTRaiseWarning(
            "[%s] The source is located in the liquid layer, "
            "therefore only the Green's Funtions for the Explosion source will be computed.\n" 
            , command);
    }

}