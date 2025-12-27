/**
 * @file   grt_lamb1.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-11
 * 
 *    求解第一类 Lamb 问题的主函数
 */

#include <stdio.h>
#include <stdlib.h>

#include "grt.h"
#include "grt/lamb/lamb1.h"
#include "grt/lamb/lamb_util.h"


/** 该子模块的参数控制结构体 */
typedef struct {
    /** 模型参数 */
    struct {
        bool active;
        real_t nu;    ///<  泊松比
    } P;

    /** 归一化时间序列 */
    struct {
        bool active;
        real_t *ts;
        int nt;
    } T;

    /** 方位角 */
    struct {
        bool active;
        real_t azimuth;  ///<  方位角，单位为度
    } A;

} GRT_MODULE_CTRL;



/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->T.ts);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/**
 * 打印使用说明
 */
static void print_help(){
printf("\n"
"[grt lamb1] %s\n\n", GRT_VERSION);printf(
"    Compute the exact closed-form solution for the first-kind Lamb problem\n"
"    (both the source and receiver are on the surface).\n"
"\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt lamb1 -P<nu> -T<t1>/<t2>/<dt> -A<azimuth>\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -P<nu>         Possion ratio of the halfspace, (0, 0.5).\n"
"\n"
"    -T<t1>/<t2>/<dt>\n"
"                   Dimensionless time.\n"
"                   <t1>: start time.\n"
"                   <t2>: end time.\n"
"                   <dt>: time interval.\n"
"\n"
"    -A<azimuth>    Azimuth in degree, from source to station.\n"
"\n"
"    -h             Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt lamb1 -P0.25 -T0/2/1e-3 -A30\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    int opt;

    while ((opt = getopt(argc, argv, ":P:T:A:h")) != -1) {
        switch (opt) {
            // 模型参数， -P<nu>
            case 'P':
                Ctrl->P.active = true;
                if(1 != sscanf(optarg, "%lf", &Ctrl->P.nu)){
                    GRTBadOptionError(P, "");
                }
                if(Ctrl->P.nu <= 0.0 || Ctrl->P.nu >= 0.5){
                    GRTBadOptionError(P, "possion ratio (%lf) is out of bound.", Ctrl->P.nu);
                }
                break;
            
            // 归一化时间序列, -Tt1/t2/dt
            case 'T':
                Ctrl->T.active = true;
                {
                    real_t a1, a2, delta;
                    if(3 != sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
                        GRTBadOptionError(T, "");
                    };
                    if(a1 < 0.0 || a2 < 0.0){
                        GRTBadOptionError(T, "t1 < 0.0 or t2 < 0.0.");
                    }
                    if(delta <= 0.0){
                        GRTBadOptionError(T, "dt <= 0.0.");
                    }
                    if(a1 > a2){
                        GRTBadOptionError(T, "t1(%f) > t2(%f).", a1, a2);
                    }

                    Ctrl->T.nt = floor((a2-a1)/delta) + 1;
                    Ctrl->T.ts = (real_t*)calloc(Ctrl->T.nt, sizeof(real_t));
                    for(int i=0; i<Ctrl->T.nt; ++i){
                        Ctrl->T.ts[i] = a1 + delta*i;
                    }
                }
                break;

            // 方位角，  -Aazimuth
            case 'A':
                Ctrl->A.active = true;
                if(1 != sscanf(optarg, "%lf", &Ctrl->A.azimuth)){
                    GRTBadOptionError(A, "");
                }
                if(Ctrl->A.azimuth < 0.0 || Ctrl->A.azimuth > 360){
                    GRTBadOptionError(A, "azimuth should be in [0, 360].");
                }
                break;
            
            GRT_Common_Options_in_Switch((char)(optopt)); 
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, P);
    GRTCheckOptionActive(Ctrl, T);
    GRTCheckOptionActive(Ctrl, A);
}



/** 模块主函数 */
int lamb1_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    // 传入参数 
    getopt_from_command(Ctrl, argc, argv);

    // 求解，输出到标准输出
    grt_solve_lamb1(Ctrl->P.nu, Ctrl->T.ts, Ctrl->T.nt, Ctrl->A.azimuth, NULL);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}