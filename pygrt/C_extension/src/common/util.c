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

#include "common/util.h"
#include "common/model.h"
#include "common/const.h"
#include "common/sacio.h"
#include "common/myfftw.h"
#include "travt/travt.h"

char ** string_split(const char *string, const char *delim, int *size)
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


const char* get_basename(const char* path) {
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
    const char *srcname, const char ch, const MYREAL delayT,
    const MYREAL wI, FFTW_HOLDER *pt_fh,
    SACHEAD *pt_hd, const char *s_output_subdir, const char *s_prefix,
    const int sgn, const MYCOMPLEX *grncplx)
{
    snprintf(pt_hd->kcmpnm, sizeof(pt_hd->kcmpnm), "%s%s%c", s_prefix, srcname, ch);
    
    char *s_outpath = NULL;
    GRT_SAFE_ASPRINTF(&s_outpath, "%s/%s.sac", s_output_subdir, pt_hd->kcmpnm);

    // 执行fft任务会修改数组，需重新置零
    reset_fftw_holder_zero(pt_fh);
    
    // 赋值复数，包括时移
    MYCOMPLEX cfac, ccoef;
    cfac = exp(I*PI2*pt_fh->df*delayT);
    ccoef = sgn;
    // 只赋值有效长度，其余均为0
    for(int i=0; i<pt_fh->nf_valid; ++i){
        pt_fh->W_f[i] = grncplx[i] * ccoef;
        ccoef *= cfac;
    }


    if(! pt_fh->naive_inv){
        // 发起fft任务 
        fftw_execute(pt_fh->plan);
    } else {
        naive_inverse_transform_double(pt_fh);
    }
    

    // 归一化，并处理虚频
    // 并转为 SAC 需要的单精度类型
    float *float_arr = (float*)malloc(sizeof(float)*pt_fh->nt);
    double fac, coef;
    coef = pt_fh->df * exp(delayT*wI);
    fac = exp(wI*pt_fh->dt);
    for(int i=0; i<pt_fh->nt; ++i){
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
 * @param[in]         pymod         模型结构体指针
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
    const char *command, const PYMODEL1D *pymod, const char *s_prefix, 
    const MYREAL wI, FFTW_HOLDER *pt_fh,
    const char *s_dist, const MYREAL dist,
    const MYREAL depsrc, const MYREAL deprcv,
    const MYREAL delayT0, const MYREAL delayV0, const bool calc_upar,
    const bool doEX, const bool doVF, const bool doHF, const bool doDC, 
    SACHEAD *pt_hd, const char *chalst,
    MYCOMPLEX *grn[SRC_M_NUM][CHANNEL_NUM], 
    MYCOMPLEX *grn_uiz[SRC_M_NUM][CHANNEL_NUM], 
    MYCOMPLEX *grn_uir[SRC_M_NUM][CHANNEL_NUM])
{
    // 文件保存子目录
    char *s_output_subdir = NULL;
    
    GRT_SAFE_ASPRINTF(&s_output_subdir, "%s_%s", s_prefix, s_dist);
    GRTCheckMakeDir(command, s_output_subdir);

    // 时间延迟 
    MYREAL delayT = delayT0;
    if(delayV0 > 0.0)   delayT += sqrt( GRT_SQUARE(dist) + GRT_SQUARE(deprcv - depsrc) ) / delayV0;
    // 修改SAC头段时间变量
    pt_hd->b = delayT;

    // 计算理论走时
    pt_hd->t0 = compute_travt1d(pymod->Thk, pymod->Va, pymod->n, pymod->isrc, pymod->ircv, dist);
    strcpy(pt_hd->kt0, "P");
    pt_hd->t1 = compute_travt1d(pymod->Thk, pymod->Vb, pymod->n, pymod->isrc, pymod->ircv, dist);
    strcpy(pt_hd->kt1, "S");

    for(int im=0; im<SRC_M_NUM; ++im){
        if(! doEX  && im==0)  continue;
        if(! doVF  && im==1)  continue;
        if(! doHF  && im==2)  continue;
        if(! doDC  && im>=3)  continue;

        int modr = SRC_M_ORDERS[im];
        int sgn=1;  // 用于反转Z分量
        for(int c=0; c<CHANNEL_NUM; ++c){
            if(modr==0 && ZRTchs[c]=='T')  continue;  // 跳过输出0阶的T分量

            // 判断是否为所需的分量
            if(strchr(chalst, ZRTchs[c]) == NULL)  continue;

            // Z分量反转
            sgn = (ZRTchs[c]=='Z') ? -1 : 1;

            write_one_to_sac(
                SRC_M_NAME_ABBR[im], ZRTchs[c], delayT, 
                wI, pt_fh,
                pt_hd, s_output_subdir, "", sgn, 
                grn[im][c]);

            if(calc_upar){
                write_one_to_sac(
                    SRC_M_NAME_ABBR[im], ZRTchs[c], delayT, 
                    wI, pt_fh,
                    pt_hd, s_output_subdir, "z", sgn*(-1), 
                    grn_uiz[im][c]);

                write_one_to_sac(
                    SRC_M_NAME_ABBR[im], ZRTchs[c], delayT, 
                    wI, pt_fh,
                    pt_hd, s_output_subdir, "r", sgn, 
                    grn_uir[im][c]);
            }

        }
    }


    GRT_SAFE_FREE_PTR(s_output_subdir);
}



void GF_freq2time_write_to_file(
    const char *command, const PYMODEL1D *pymod, 
    const char *s_output_dir, const char *s_modelname, const char *s_depsrc, const char *s_deprcv,    
    const MYREAL wI, FFTW_HOLDER *pt_fh,
    const MYINT nr, char *s_dists[nr], const MYREAL dists[nr], MYREAL travtPS[nr][2],
    const MYREAL depsrc, const MYREAL deprcv,
    const MYREAL delayT0, const MYREAL delayV0, const bool calc_upar,
    const bool doEX, const bool doVF, const bool doHF, const bool doDC, 
    const char *chalst,
    MYCOMPLEX *grn[nr][SRC_M_NUM][CHANNEL_NUM], 
    MYCOMPLEX *grn_uiz[nr][SRC_M_NUM][CHANNEL_NUM], 
    MYCOMPLEX *grn_uir[nr][SRC_M_NUM][CHANNEL_NUM])
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
    hd.user1 = pymod->Va[pymod->ircv];
    hd.user2 = pymod->Vb[pymod->ircv];
    hd.user3 = pymod->Rho[pymod->ircv];
    hd.user4 = RONE/pymod->Qa[pymod->ircv];
    hd.user5 = RONE/pymod->Qb[pymod->ircv];
    // 写入震源点的Vp,Vs,rho
    hd.user6 = pymod->Va[pymod->isrc];
    hd.user7 = pymod->Vb[pymod->isrc];
    hd.user8 = pymod->Rho[pymod->isrc];

    char *s_output_dirprefx = NULL;
    GRT_SAFE_ASPRINTF(&s_output_dirprefx, "%s/%s_%s_%s", s_output_dir, s_modelname, s_depsrc, s_deprcv);
    
    // 做反傅里叶变换，保存SAC文件
    for(int ir=0; ir<nr; ++ir){
        hd.dist = dists[ir];

        single_freq2time_write_to_file(
            command, pymod, s_output_dirprefx, 
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
    if(pymod->Vb[pymod->isrc]==0.0){
        fprintf(stderr, "[%s] " BOLD_YELLOW 
            "The source is located in the liquid layer, "
            "therefore only the Green's Funtions for the Explosion source will be computed.\n" 
            DEFAULT_RESTORE, command);
    }

}