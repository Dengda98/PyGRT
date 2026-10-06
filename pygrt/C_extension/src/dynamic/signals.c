/**
 * @file   signals.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12-2
 * 
 *    时间函数生成与信号处理
 *    时间函数及信号处理计算使用 real_t，存储精度由 I/O 层处理
 * 
 *    信号长度应能整除采样间隔。
 * 
 *    时间函数按矩形法进行面积归一化，以匹配离散卷积（雷克子波保留单位峰值）
 *    自定义时间函数检查离散积分，非单位面积时警告并自动归一化
 * 
 */

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "grt/dynamic/signals.h"
#include "grt/common/util.h"

#include "grt/common/checkerror.h"

/**
 * 按矩形法将时间函数的面积归一化，返回归一化前的面积
 *
 * @param[in,out]  x   时间函数样本
 * @param[in]      nx  样本数
 * @param[in]      dt  采样间隔
 */
static real_t grt_normalize_time_function(real_t *x, int nx, real_t dt)
{
    // 使用与离散卷积相同的 dt 乘样本和，并检查正负幅值抵消
    real_t area = 0.0, absarea = 0.0;
    for(int n=0; n<nx; ++n) {
        area    += x[n] * dt;
        absarea += fabs(x[n]) * dt;
    }
    if(absarea == 0.0 || GRT_ISCLOSE(area / absarea, 0.0)) {
        GRTRaiseError("Time function has zero or numerically singular integral.");
    }
    for(int n=0; n<nx; ++n) {
        x[n] /= area;
    }
    return area;
}

/**
 * 检查不含破裂延迟的时间函数参数
 * @param[in] tftype    时间函数类型
 * @param[in] tfparams  时间函数参数或自定义文件路径
 */
static bool check_time_function_base(const char tftype, const char *tfparams){

    // 脉冲
    if(GRT_SIG_IMPULSE == tftype) {
        return !*tfparams;
    }
    // 抛物波
    else if(GRT_SIG_PARABOLA == tftype){
        real_t t0 = 0.0;
        if(1 != sscanf(tfparams, "%lf", &t0)) {
            return false;
        }
        if(t0 <= 0) {
            GRTRaiseError("t0(%s) should be larger than 0.\n", tfparams);
        }
    }
    // 梯形波
    else if(GRT_SIG_TRAPEZOID == tftype){
        real_t t1 = 0.0, t2 = 0.0, t3 = 0.0;
        if(3 != sscanf(tfparams, "%lf/%lf/%lf", &t1, &t2, &t3)) {
            return false;
        }
        if(t1 < 0.0 || t2 < 0.0 || t3 <= 0.0){
            GRTRaiseError("It should be t1>=0.0, t2>=0.0 and t3>0.0 (%s).\n", tfparams);
        }
        if(t1 > t2 || t2 > t3) {
            GRTRaiseError("It should be t1<=t2<=t3 (%s).\n", tfparams);
        }
    }
    // 雷克子波
    else if(GRT_SIG_RICKER == tftype){
        real_t f0;
        if(1 != sscanf(tfparams, "%lf", &f0)) {
            return false;
        }
        if(f0 <= 0) {
            GRTRaiseError("f0(%s) should be larger than 0.\n", tfparams);
        }
    }
    // 自定义时间函数
    else if(GRT_SIG_CUSTOM == tftype){
        // tfparams 为自定义时间函数的文件名，检查文件是否存在
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

/**
 * 分离时间函数参数与非负破裂延迟，返回独立参数副本
 * @param[in]  params  含可选 +d 后缀的时间函数参数
 * @param[out] delay   非负破裂延迟，s
 */
static char *split_time_delay(const char *params, real_t *delay)
{
    char *base = strdup(params ? params : "");
    *delay = 0.0;

    char *suffix = strstr(base, "+d");
    if(suffix) {
        char extra;
        if(sscanf(suffix + 2, "%lf%c", delay, &extra) != 1 || *delay < 0.0) {
            GRTRaiseError("Time-function delay must be nonnegative.");
        }
        *suffix = '\0';
    }
    return base;
}

bool grt_check_tftype_tfparams(const char tftype, const char *tfparams)
{
    real_t delay;
    char *base = split_time_delay(tfparams, &delay);
    bool good = check_time_function_base(tftype, base);
    free(base);
    return good;
}

real_t * grt_get_time_function(int *TFnt, real_t dt, const char tftype, const char *tfparams)
{
    if(dt <= 0.0) {
        GRTRaiseError("Invalid time-function sampling interval.");
    }

    real_t delay;
    char *base = split_time_delay(tfparams, &delay);
    if(!check_time_function_base(tftype, base)) {
        GRTRaiseError("Invalid time function.");
    }

    // 按类型生成时间函数
    real_t *tfarr = NULL;
    int tfnt = 0;
    if(GRT_SIG_IMPULSE == tftype){
        tfnt = 1;
        tfarr = GRT_SAFE_CALLOC(tfnt, sizeof(*tfarr));
        tfarr[0] = 1.0 / dt;
    }
    // 抛物波
    else if(GRT_SIG_PARABOLA == tftype){
        real_t t0 = 0.0;
        sscanf(base, "%lf", &t0);
        tfarr = grt_get_parabola_wave(dt, &t0, &tfnt);
    }

    // 梯形波
    else if(GRT_SIG_TRAPEZOID == tftype){
        real_t t1 = 0.0, t2 = 0.0, t3 = 0.0;
        sscanf(base, "%lf/%lf/%lf", &t1, &t2, &t3);
        tfarr = grt_get_trap_wave(dt, &t1, &t2, &t3, &tfnt);
    }

    // 雷克子波
    else if(GRT_SIG_RICKER == tftype){
        real_t f0 = 0.0;
        sscanf(base, "%lf", &f0);
        tfarr = grt_get_ricker_wave(dt, f0, &tfnt);
    }

    // 自定义时间函数
    else if(GRT_SIG_CUSTOM == tftype){
        tfarr = grt_get_custom_wave(&tfnt, base);

        // 自定义时间函数共用矩形法归一化，原始面积不为 1 时提示用户
        real_t area = grt_normalize_time_function(tfarr, tfnt, dt);
        if(!GRT_ISCLOSE(area, 1.0)) {
            GRTRaiseWarning("Custom time function sequence sum is %.7g, expected %.7g (1/dt); normalizing automatically.",
                area/dt, 1.0/dt);
        }
    }

    *TFnt = tfnt;
    free(base);

    // 延迟取最近采样点，半采样点向上取整，再通过前导零平移时间函数
    int shift = (int)llround(delay / dt);
    if(shift) {
        real_t *padded = GRT_SAFE_CALLOC(*TFnt + shift, sizeof(*padded));
        memcpy(padded + shift, tfarr, *TFnt * sizeof(*tfarr));
        free(tfarr);
        tfarr = padded;
        *TFnt += shift;
    }
    return tfarr;
}

real_t *grt_time_function_from_option(const char *option, real_t dt, int *nt)
{
    // NULL 使用脉冲，否则解析完整 -D 选项
    if(!option) {
        option = "-Di";
    }
    if(strncmp(option, "-D", 2) != 0 || !option[2]) {
        GRTRaiseError("Expected a complete -D time-function option.");
    }

    char type = option[2];
    const char *params = option + 3;
    if(*params == '/') {
        ++params;
        if(!*params || strncmp(params, "+d", 2) == 0) {
            GRTRaiseError("Time-function parameters must follow '/'.");
        }
    } else if(*params && strncmp(params, "+d", 2) != 0) {
        GRTRaiseError("Expected -Dtftype[/tfparams][+d<delay>].");
    }

    return grt_get_time_function(nt, dt, type, params);
}



void grt_oaconvolve(const real_t *x, int nx, const real_t *h, int nh, real_t *y, int ny, bool iscircular) {
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


real_t grt_trap_area(const real_t *x, int nx, real_t dt){
    real_t area = 0.0;
    for(int i=0; i<nx-1; ++i){
        area += (x[i] + x[i+1]) * 0.5 * dt;
    }
    return area;
}


void grt_trap_integral(real_t *x, int nx, real_t dt){
    // 矩形法
    // x[0] = 0.0; // 边界条件
    // for(int i=1; i<nx; ++i){
    //     x[i] = x[i]*dt + x[i-1];
    // }
    // 梯形法
    real_t lastx=x[0], tmp;
    x[0] = 0.0;
    for(int i=1; i<nx; ++i){
        tmp = x[i];
        x[i] = 0.5*(x[i] + lastx)*dt + x[i-1];
        lastx = tmp;
    }
}



void grt_differential(real_t *x, int nx, real_t dt){
    if(nx == 1) {
        x[0] = 0.0;
        return;
    }
    // 中心差分
    real_t tmp, x0=x[0];
    real_t h=2.0*dt;
    x[0] = (x[1]-x0)/dt;
    for(int i=1; i<nx-1; ++i){
        tmp = (x[i+1] - x0)/h;
        x0 = x[i];
        x[i] = tmp;
    }
    x[nx-1] = (x[nx-1] - x0)/dt;
}


real_t * grt_get_parabola_wave(real_t dt, real_t *Tlen, int *Nt){
    if(dt <= 0.0 || *Tlen <= 0.0) {
        GRTRaiseError("Parabolic duration and sampling interval must be positive.");
    }

    // 截止时刻向上对齐到采样网格，至少保留一个位于两个零端点之间的样本
    real_t sample = *Tlen / dt;
    real_t nearest = round(sample);
    int last = GRT_ISCLOSE(*Tlen, nearest * dt) ? (int)nearest : (int)ceil(sample);
    if(last < 2) {
        GRTRaiseError("Window length of time function is too short.");
    }

    // 归一化时间从 0 到 1，对应抛物线的两个零端点
    int nt = last + 1;
    real_t *arr = GRT_SAFE_CALLOC(nt, sizeof(*arr));
    for(int n=0; n<nt; ++n) {
        real_t phase = (real_t)n / last;
        arr[n] = 4.0 * phase * (1.0 - phase);
    }

    grt_normalize_time_function(arr, nt, dt);

    *Tlen = last * dt;
    *Nt = nt;
    return arr;
}


real_t * grt_get_trap_wave(real_t dt, real_t *T1, real_t *T2, real_t *T3, int *Nt){
    real_t times[3] = {*T1, *T2, *T3};
    if(dt <= 0.0) {
        GRTRaiseError("Invalid time-function sampling interval.");
    }
    if(times[0] < 0.0 || times[0] > times[1] || times[1] > times[2] || times[2] <= 0.0) {
        GRTRaiseError("Trapezoidal cutoffs must satisfy 0 <= t1 <= t2 <= t3 and t3 > 0.");
    }

    // 三个截止时刻向上对齐到采样网格，已接近网格点的时刻只消除浮点误差
    int indices[3];
    for(int k=0; k<3; ++k) {
        real_t sample = times[k] / dt;
        real_t nearest = round(sample);
        indices[k] = GRT_ISCLOSE(times[k], nearest * dt) ? (int)nearest : (int)ceil(sample);
    }

    // i1、i2、i3 分别为上坡、平台、下坡的截止下标，零时长段共享截止下标
    int i1 = *T1 > 0.0 ? GRT_MAX(1, indices[0]) : 0;
    int i2 = GRT_ISCLOSE(*T1, *T2) ? i1 : GRT_MAX(i1, indices[1]);
    int i3 = GRT_ISCLOSE(*T2, *T3) ? i2 : GRT_MAX(i2 + 1, indices[2]);

    // 时间函数至少保留一个采样间隔
    if(i3 == 0) {
        i2 = i3 = 1;
    }

    int nt = i3 + 1;
    real_t *arr = GRT_SAFE_CALLOC(nt, sizeof(*arr));
    for(int n=0; n<nt; ++n) {
        if(n < i1) {
            arr[n] = (real_t)n / i1;                   // 上坡：从 0 线性升至 1
        } else if(n <= i2) {
            arr[n] = 1.0;                              // 平台：三角波在这里仅有一个顶点
        } else {
            arr[n] = (real_t)(i3 - n) / (i3 - i2);      // 下坡：从 1 线性降至 0
        }
    }

    grt_normalize_time_function(arr, nt, dt);

    *T1 = i1 * dt;
    *T2 = i2 * dt;
    *T3 = i3 * dt;
    *Nt = nt;
    return arr;
}


real_t * grt_get_ricker_wave(real_t dt, real_t f0, int *Nt){
    if(dt <= 0.0 || f0 <= 0.0) {
        GRTRaiseError("Ricker frequency and sampling interval must be positive.");
    }
    if(1.0 / dt <= 2.0 * f0) {
        GRTRaiseError("Compare to sampling freq (%.3f), dominant freq (%.3f) is too high.", 1.0 / dt, f0);
    }

    // 峰值时刻为 1/f0，采样窗口覆盖其前后各一个周期
    real_t peak_time = 1.0 / f0;
    int nt = 2 * ((int)floor(peak_time / dt) + 1);
    real_t *arr = GRT_SAFE_CALLOC(nt, sizeof(*arr));

    // 按下标直接计算时间，避免累加 dt 的舍入误差，保留解析波形的单位峰值
    for(int n=0; n<nt; ++n) {
        real_t phase = PI * f0 * ((real_t)n * dt - peak_time);
        real_t phase2 = phase * phase;
        arr[n] = (1.0 - 2.0 * phase2) * exp(-phase2);
    }

    *Nt = nt;
    return arr;
}


real_t * grt_get_custom_wave(int *Nt, const char *tfparams){
    FILE *fp = fopen(tfparams, "r");
    if(fp == NULL) {
        GRTRaiseError("Custom time function file open error.");
    }

    // 逐行读取一列振幅，跳过空行和注释行
    real_t *tfarr = NULL;
    char *line = NULL;
    size_t len = 0;
    size_t lineno = 0;
    int nt = 0;
    while(grt_getline(&line, &len, fp) != -1) {
        lineno++;
        if(grt_is_comment_or_empty_line(line)) {
            continue;
        }

        // 每个非注释行只能包含一列振幅值
        real_t value = 0.0;
        char extra = '\0';
        if(sscanf(line, " %lf %c", &value, &extra) != 1) {
            GRTRaiseError("custom time function file should contain exactly one column at line %zu.\n", lineno);
        }

        tfarr = GRT_SAFE_REALLOC(tfarr, (nt + 1) * sizeof(*tfarr));
        tfarr[nt] = value;
        nt++;
    }

    if(nt == 0) {
        GRTRaiseError("Custom time function file is empty.");
    }

    fclose(fp);
    GRT_SAFE_FREE_PTR(line);

    *Nt = nt;
    return tfarr;
}


void grt_free1d(void *pt){
    GRT_SAFE_FREE_PTR(pt);
}
