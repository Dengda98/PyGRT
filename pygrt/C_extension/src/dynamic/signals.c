/**
 * @file   signals.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12-2
 * 
 *    常见的时间函数，为了sac的兼容性，以下均使用float指针，
 *    并在函数内部使用malloc申请内存。
 * 
 *    信号长度应能整除采样间隔。
 * 
 *    所有时间函数的最大值都是1.0
 * 
 */

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "grt/dynamic/signals.h"
#include "grt/common/const.h"
#include "grt/common/util.h"

#include "grt/common/checkerror.h"


bool grt_check_tftype_tfparams(const char tftype, const char *tfparams){

    // 抛物波
    if(GRT_SIG_PARABOLA == tftype){
        float t0=0.0;
        if(1 != sscanf(tfparams, "%f", &t0))  return false;
        if(t0 <= 0){
            GRTRaiseError("t0(%s) should be larger than 0.\n", tfparams);
        }
    }
    // 梯形波
    else if(GRT_SIG_TRAPEZOID == tftype){
        float t1=0.0, t2=0.0, t3=0.0;
        if(3 != sscanf(tfparams, "%f/%f/%f", &t1, &t2, &t3))   return false;
        if(t1 < 0.0 || t2 < 0.0 || t3 <= 0.0){
            GRTRaiseError("It should be t1>=0.0, t2>=0.0 and t3>0.0 (%s).\n", tfparams);
        }
        if(! (t1 <= t2 && t2 < t3)){
            GRTRaiseError("It should be t1<=t2<t3 (%s).\n", tfparams);
        }
    }
    // 雷克子波
    else if(GRT_SIG_RICKER == tftype){
        float f0;
        if(1 != sscanf(tfparams, "%f", &f0))  return false;
        if(f0 <= 0){
            GRTRaiseError("f0(%s) should be larger than 0.\n", tfparams);
        }
    }
    // 自定义时间函数
    else if(GRT_SIG_CUSTOM == tftype){
        // tfparams为存储自定义时间函数的文件名
        // 检查文件是否存在
        if(access(tfparams, F_OK) != 0){
            GRTRaiseError("(%s) not exists.\n", tfparams);
        }
    }
    // 不符合要求
    else{
        GRTRaiseError("Unsupported time function type '%c'.\n", tftype);
    }

    return true;
}


float * grt_get_time_function(int *TFnt, float dt, const char tftype, const char *tfparams){
    // 获得时间函数
    float *tfarr=NULL;
    int tfnt=0;
    // 抛物波
    if(GRT_SIG_PARABOLA == tftype){
        float t0=0.0;
        sscanf(tfparams, "%f", &t0);
        tfarr = grt_get_parabola_wave(dt, &t0, &tfnt);
    }
    // 梯形波
    else if(GRT_SIG_TRAPEZOID == tftype){
        float t1=0.0, t2=0.0, t3=0.0;
        sscanf(tfparams, "%f/%f/%f", &t1, &t2, &t3);
        tfarr = grt_get_trap_wave(dt, &t1, &t2, &t3, &tfnt);
    }
    // 雷克子波
    else if(GRT_SIG_RICKER == tftype){
        float f0=0.0;
        sscanf(tfparams, "%f", &f0);
        tfarr = grt_get_ricker_wave(dt, f0, &tfnt);
    }
    // 自定义时间函数
    else if(GRT_SIG_CUSTOM == tftype){
        tfarr = grt_get_custom_wave(&tfnt, tfparams);
    }

    *TFnt = tfnt;
    return tfarr;
}



void linear_convolve_time_function(float *arr, int nt, float dt, const char tftype, const char *tfparams, float **TFarr, int *TFnt){
    // 获得时间函数
    float *tfarr=NULL;
    int tfnt=0;
    tfarr = grt_get_time_function(&tfnt, dt, tftype, tfparams);

    float *yarr = (float*)calloc(nt, sizeof(float));
    // 线性卷积
    grt_oaconvolve(arr, nt, tfarr, tfnt, yarr, nt, false);

    // 原地更改
    for(int i=0; i<nt; ++i){
        arr[i] = yarr[i] * dt; // dt为卷积的系数
    }
    GRT_SAFE_FREE_PTR(yarr);


    if(TFarr!=NULL || TFnt!=NULL){
        if(TFarr!=NULL) *TFarr = tfarr;
        if(TFnt!=NULL)  *TFnt  = tfnt;
    } else {
        GRT_SAFE_FREE_PTR(tfarr);
    }   

}


void grt_oaconvolve(float *x, int nx, float *h, int nh, float *y, int ny, bool iscircular) {
    if(iscircular){
        for(int n=0; n<ny; ++n) {
            y[n] = 0.0;
            for(int k=0; k<nh; ++k) {
                y[n] += x[(n - k + nx)%nx] * h[k];
            }
        }
    } else {
        for(int n=0; n<ny; ++n) {
            y[n] = 0.0;
            for(int k=0; k<nh; ++k) {
                if (n - k >= 0 && n - k < nx) {
                    y[n] += x[n - k] * h[k]; // 计算卷积值
                }
            }
        }
    }
}


float grt_trap_area(const float *x, int nx, float dt){
    float area = 0.0;
    for(int i=0; i<nx-1; ++i){
        area += (x[i] + x[i+1])*0.5/dt;
    }
    return area;
}


void grt_trap_integral(float *x, int nx, float dt){
    // 矩形法
    // x[0] = 0.0; // 边界条件
    // for(int i=1; i<nx; ++i){
    //     x[i] = x[i]*dt + x[i-1];
    // }
    // 梯形法
    float lastx=x[0], tmp;
    x[0] = 0.0;
    for(int i=1; i<nx; ++i){
        tmp = x[i];
        x[i] = 0.5*(x[i] + lastx)*dt + x[i-1];
        lastx = tmp;
    }
}



void grt_differential(float *x, int nx, float dt){
    // 中心差分
    float tmp, x0=x[0];
    float h=2.0*dt;
    x[0] = (x[1]-x0)/dt;
    for(int i=1; i<nx-1; ++i){
        tmp = (x[i+1] - x0)/h;
        x0 = x[i];
        x[i] = tmp;
    }
    x[nx-1] = (x[nx-1] - x0)/dt;
}


float * grt_get_parabola_wave(float dt, float *Tlen, int *Nt){
    float tlen = *Tlen;
    int nt = floorf(tlen/dt);
    if(fabsf(tlen - nt*dt) <= 1e-6) nt--;
    if(nt==0) {
        GRTRaiseError("window length of time function is too short.\n");
    }
    nt += 2;
    tlen = (nt-1)*dt;

    float *arr = (float*)calloc(nt, sizeof(float));
    float fac=1.0/(nt-1);
    float pha=0.0;
    for(int n=0; n<nt; ++n){
        arr[n] = 1.0 - (pha - 0.5)*(pha - 0.5)*4.0;
        pha += fac;
    }

    *Tlen = tlen;
    *Nt = nt;
    return arr;
}


float * grt_get_trap_wave(float dt, float *T1, float *T2, float *T3, int *Nt){
    // 如果t1==t2，则退化为三角波
    bool istriangle = (*T1==*T2)? true : false;

    // 微调T1
    float t1 = *T1;
    int n1 = floorf(t1/dt);
    if((fabsf(t1 - n1*dt)) <= 1e-6) n1--;
    n1 += 2;
    t1 = (n1-1)*dt;
    // 微调T2
    float t2;
    int n2;
    if(istriangle){
        n2 = 1;
        t2 = t1;
    } else {
        t2 = *T2;
        n2 = floorf((t2-t1)/dt);
        if((fabsf(t2-t1 - n2*dt)) <= 1e-6) n2--;
        n2 += 2;
        t2 = t1 + (n2-1)*dt;
    }
    
    // 微调T3
    float t3 = *T3;
    int n3 = floorf((t3-t2)/dt);
    if((fabsf(t3-t2 - n3*dt)) <= 1e-6) n3--;
    n3 += 2;
    t3 = t2 + (n3-1)*dt;

    // 总点数
    int nt = n1+n2+n3 - 2;
    float *arr = (float*)calloc(nt, sizeof(float));

    float fac=0.0, y=0.0;
    // 上坡
    fac = 1.0/(n1-1);
    y = 0.0;
    for(int n=0; n<n1; ++n){
        arr[n] = y;
        y += fac;
    }

    // 平台
    for(int n=n1-1; n<n1+n2-1; ++n){
        arr[n] = 1.0;
    }

    // 下坡
    fac = 1.0/(n3-1);
    y = 1.0;
    for(int n=n1+n2-2; n<n1+n2+n3-2; ++n){
        arr[n] = y;
        y -= fac;
    }
    

    *T1 = t1;
    *T2 = t2;
    *T3 = t3;
    *Nt = nt;
    return arr;
}


float * grt_get_ricker_wave(float dt, float f0, int *Nt){
    if(1.0/dt <= 2.0*f0) { // 在当前采样率下，主频f0过高
        GRTRaiseError("Compare to sampling freq (%.3f), dominant freq (%.3f) is too high.\n", 1.0/dt, f0);
    }

    float t0 = 1.0/f0;
    int nt = (floorf(t0/dt) + 1) * 2; // 估计2倍长度够包含
    float *arr = (float*)calloc(nt, sizeof(float));

    float PPI = PI*PI;
    float ff0 = f0*f0;
    float a, t=0.0;
    for(int i=0; i<nt; ++i){
        a = PPI*ff0*(t-t0)*(t-t0);
        arr[i] = (1.0 - 2.0*a)*expf(-a);
        t += dt;
    }

    *Nt = nt;
    return arr;
}


float * grt_get_custom_wave(int *Nt, const char *tfparams){
    float *tfarr = (float*)malloc(sizeof(float)*1);
    FILE *fp;
    if((fp = fopen(tfparams, "r")) == NULL){
        GRTRaiseError("custom time function file open error.\n");
    }

    // 逐行读入
    size_t len;
    char *line = NULL;

    int nt = 0;
    while(grt_getline(&line, &len, fp) != -1) {
        // 注释行
        if(grt_is_comment_or_empty(line))  continue;

        tfarr = (float*)realloc(tfarr, sizeof(float)*(nt+1));
        if(sscanf(line, " %f", &tfarr[nt]) < 1){
            GRTRaiseError("custom time function file read error.\n");
        }
        nt++;
    }

    if(nt == 0){
        GRTRaiseError("custom time function file read error. Empty?\n");
    }

    fclose(fp);
    GRT_SAFE_FREE_PTR(line);

    *Nt = nt;
    return tfarr;
}


void grt_free1d(void *pt){
    GRT_SAFE_FREE_PTR(pt);
}