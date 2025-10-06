/**
 * @file   grt_travt.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-08
 * 
 *    计算一维均匀半无限层状介质的初至走时，思路借鉴了CPS330程序包的time96.f
 * 
 */

#include "grt/travt/travt.h"
#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/common/util.h"

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    char *name;
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;
        GRT_MODEL1D *mod1d;         ///< 模型结构体指针
    } M;
    /** 震源和接收器深度 */
    struct {
        bool active;
        MYREAL depsrc;
        MYREAL deprcv;
        char *s_depsrc;
        char *s_deprcv;
    } D;
    /** 震中距 */
    struct {
        bool active;
        char **s_rs;
        MYREAL *rs;
        MYINT nr;
    } R;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->name);
    GRT_SAFE_FREE_PTR(Ctrl->M.s_modelpath);
    grt_free_mod1d(Ctrl->M.mod1d);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_depsrc);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_deprcv);

    GRT_SAFE_FREE_PTR(Ctrl->R.rs);
    GRT_SAFE_FREE_PTR_ARRAY(Ctrl->R.s_rs, Ctrl->R.nr);

    GRT_SAFE_FREE_PTR(Ctrl);
}



MYREAL grt_compute_travt1d(
    const MYREAL *Thk, const MYREAL *Vel0, const int nlay, 
    const int isrc, const int ircv, const MYREAL dist)
{
    // 以防速度数组中存在零速度的情况，这里新建数组以去除0速度
    MYREAL *Vel = (MYREAL*)malloc(sizeof(MYREAL)*nlay);
    for(int i=0; i<nlay; ++i){
        Vel[i] = (Vel0[i] <= 0.0)? 1e-6 : Vel0[i];  // 给一个极慢值
    }

    // 根据互易原则，震源和场点可交换
    int imin, imax;
    imin = (isrc < ircv)? isrc : ircv;
    imax = (isrc < ircv)? ircv : isrc;

    // 如果震源与台站之间存在零速度层，则此时单一震相类型的走时不必再计算，
    // 直接返回默认值
    for(int i = imin; i <= imax; ++i){
        if(Vel0[i] == 0.0){
            GRT_SAFE_FREE_PTR(Vel);
            return -12345.00;
        }
    }
    

    // 震源和场点速度
    MYREAL vsrc = Vel[isrc];
    MYREAL vrcv = Vel[ircv];
    MYREAL vmax = (vsrc < vrcv)? vrcv : vsrc;
    MYREAL vmin = (vsrc < vrcv)? vsrc : vrcv;
    // 震源和场点深度
    MYREAL depsrc=0.0, deprcv=0.0;
    for(int i=0; i<isrc; ++i) {depsrc += Thk[i];} 
    for(int i=0; i<ircv; ++i) {deprcv += Thk[i];} 
    MYREAL depdif = fabs(depsrc - deprcv);


    // 初始化走时
    MYREAL travt = 9.0e30;

    //=====================================================    
    // 对于四种情况逐一讨论，取最小走时
    // + 直达波
    //    - 两点位于同一层位，直接用距离除速度
    //    - 两点位于不同层位，使用二分迭代找到最佳慢度
    // + 透射波
    //    - 射线从震源向上出发
    //    - 射线从震源向下出发
    //    射线会包括一个单边和对称的双边。
    // 
    //=====================================================    

    /*
    //--------------X----------------------
    //  imin         \ 
    //-------------------------------------
    //                 \ 
    //-------------------------------------
    //                   ...
    //-------------------------------------
    //                      \ 
    //-----------------------X-------------
    //   imax
    //-------------------------------------
    //               ...
    //-------------------------------------*/
    //           (halfspace)
    //
    // =========================================================
    // ------------------- 同层直达波 ----------------------
    if(abs(isrc - ircv)==1){ // 位于同一物理层
        travt = sqrt(dist*dist + depdif*depdif) / vsrc;
        // printf("direct wave in same layer, travt=%f\n", travt);
    }
    else {
    // ------------------- 不同层直达波 ----------------------
    // -------------- 使用二分迭代法进行打靶 ------------------
        // 最大迭代次数
        const int nloop=50; 
        // 最小震中距差
        const double minX=1e-3;
        // 找到慢度上限，准确说是各层中最大慢度的最小值
        MYREAL pmax0=1.0/vmax;
        for(int i=imin; i<=imax; ++i){
            if(Thk[i] == 0.0) continue;
            if(pmax0 > 1.0/Vel[i])  pmax0 = 1.0/Vel[i];
        }
        // 初始化一些迭代变量
        MYREAL pmin=0.0, pmax=pmax0;
        MYREAL p;
        MYREAL s, c, v, h;
        MYREAL x = 0.0;
        MYREAL t = 0.0;
        MYREAL tint = 0.0;
        MYREAL dxdp = 0.0;
        for(int iter=0; iter<nloop; ++iter){
            x = t = tint = dxdp = 0.0;
            p = (pmin+pmax)/2.0;
            for(int i=imin; i<imax; ++i){
                h = Thk[i];
                if(h == 0.0) continue;
                v = Vel[i];
                s = p*v;
                c = sqrt(1.0 - s*s);
                t += h/(v*c);
                x += h*s/c;
                dxdp += h*v/(c*c*c);
                // printf("i=%d, t=%f, x=%f\n", i, t, x);
            }

            if(x < dist){
                pmin = p;
            } else if(x > dist){
                pmax = p;
            } else {
                break;
            }

            if(fabs(x - dist) < minX) break;


            // printf("iter=%d, t=%f\n", iter, t);
        } // 结束迭代
        
        travt = t;

        // printf("direct wave in different layer, travt=%f\n", travt);
    }

    /*
    //---------------------------------------------------------------
    //                             ...
    //---------------------------------------------------------------
    //                         ____..._____    
    //---------------------------------------------------------------
    //                       /              \ 
    //---------------------------------------------------------------
    //                   ...                 ...  
    //---------------------------------------------------------------
    //                   /                     \    
    //------------------------------------------X-------------------
    //    imin         /
    //---------------------------------------------------------------
    //              ...
    //---------------------------------------------------------------
    //             /
    //------------X--------------------------------------------------
    //    imax
    //---------------------------------------------------------------
    //                           ...
    //---------------------------------------------------------------*/
    //                        (halfspace)
    //
    // 
    //=====================================================================
    //------------------- 向上出射的射线，考虑透射 -----------------
    if(Thk[0] > 0.0){  // 存在射线向上的基本条件
        MYREAL v, p, h, c;
        MYREAL sumt, sumx;
        bool badrefrac = false;
        // 找到透射位置
        for(int m=imin-1; m>=0; --m){
            h = Thk[m];
            if(h == 0.0) continue;
            v = Vel[m];
            p = 1.0/v;
            badrefrac = false;

            // 两点处的速度必须比透射点的速度低，且透射点速度是整个路径上最快速度
            if(vmin >= v)  continue;
            if(vmax >= v)  continue;

            sumt = sumx = 0.0;
            // imax到imin的单边
            for(int i=imin; i<imax; ++i){
                if(Vel[i] > v) {
                    badrefrac = true;
                    break;
                } 
                c = sqrt(fabs(1.0 - p*p*Vel[i]*Vel[i]));
                sumt += Thk[i]/(Vel[i]*c);
                // 走时已经超过目前最小走时，不必再讨论这一层的透射
                if(sumt > travt){
                    badrefrac = true;
                    break;
                }

                sumx += Thk[i]*p*Vel[i]/c;
                // 理论震中距已超过，不必再讨论这一层的透射
                if(sumx > dist){
                    badrefrac = true;
                    break;
                }
            }

            // 不考虑透射部分，走时已经超过当前最小走时，透射循环可结束
            if(sumt > travt) break;

            if(badrefrac) continue;

            // m到imin的双边
            for(int i=m+1; i<imin; ++i){
                if(Vel[i] > v) {
                    badrefrac = true;
                    break;
                } 
                c = sqrt(fabs(1.0 - p*p*Vel[i]*Vel[i]));
                sumt += 2.0*Thk[i]/(Vel[i]*c);
                // 走时已经超过目前最小走时，不必再讨论这一层的透射
                if(sumt > travt){
                    badrefrac = true;
                    break;
                }

                // 理论震中距已超过，不必再讨论这一层的透射
                sumx += 2.0*Thk[i]*p*Vel[i]/c;
                if(sumx > dist){
                    badrefrac = true;
                    break;
                }
            }

            // 不考虑透射部分，走时已经超过当前最小走时，透射循环可结束
            if(sumt > travt) break;

            if(badrefrac) continue;

            // printf("up m=%d, refracted wave, sumt=%f, sumx=%f\n", m, sumt, sumx);

            // 统计走时 
            if(dist >= sumx){
                sumt += (dist - sumx)/v;
                if(sumt < travt){
                    travt = sumt;
                    // printf("refracted wave in layer %d, travt=%f\n", m, travt);
                }
            }

        }  // END 寻找投射位置

    } // END 射线向上传的讨论


    /*
    //-------------------------------------------------------------
    //                             ...
    //-------------------------------------------------------------
    //
    //-------X-----------------------------------------------------
    //  imin  \ 
    //-------------------------------------------------------------
    //          ...
    //-------------------------------------------------------------
    //            \ 
    //----------------------------------------X--------------------
    //  imax        \                        /
    //-------------------------------------------------------------
    //               ...                 ...
    //-------------------------------------------------------------
    //                  \                 /
    //-------------------------------------------------------------
    //                    ‾‾‾‾‾‾ ... ‾‾‾‾
    //-------------------------------------------------------------
    //                           ...
    //---------------------------------------------------------------*/
    //                        (halfspace)
    //
    //===================================================================
    //------------------- 向下出射的射线，考虑透射 ----------------- 
    // 找到透射位置
    for(int m=imax+1; m<nlay; ++m){
        MYREAL v, p, h, c;
        MYREAL sumt, sumx;
        bool badrefrac = false;
        h = Thk[m];
        if(h == 0.0) continue;
        v = Vel[m];
        p = 1.0/v;
        badrefrac = false;

        // 两点处的速度必须比透射点的速度低，且透射点速度是整个路径上最快速度
        if(vmin >= v)  continue;
        if(vmax >= v)  continue;

        sumt = sumx = 0.0;
        // imax到imin的单边
        for(int i=imin; i<imax; ++i){
            if(Vel[i] > v) {
                badrefrac = true;
                break;
            } 
            c = sqrt(fabs(1.0 - p*p*Vel[i]*Vel[i]));
            sumt += Thk[i]/(Vel[i]*c);
            // 走时已经超过目前最小走时，不必再讨论这一层的透射
            if(sumt > travt){
                badrefrac = true;
                break;
            }

            sumx += Thk[i]*p*Vel[i]/c;
            // 理论震中距已超过，不必再讨论这一层的透射
            if(sumx > dist){
                badrefrac = true;
                break;
            }
        }

        // 不考虑透射部分，走时已经超过当前最小走时，透射循环可结束
        if(sumt > travt) break;

        if(badrefrac) continue;

        // m到imin的双边
        for(int i=imax; i<m; ++i){
            if(Vel[i] > v) {
                badrefrac = true;
                break;
            } 
            c = sqrt(fabs(1.0 - p*p*Vel[i]*Vel[i]));
            sumt += 2.0*Thk[i]/(Vel[i]*c);
            // 走时已经超过目前最小走时，不必再讨论这一层的透射
            if(sumt > travt){
                badrefrac = true;
                break;
            }

            // 理论震中距已超过，不必再讨论这一层的透射
            sumx += 2.0*Thk[i]*p*Vel[i]/c;
            if(sumx > dist){
                badrefrac = true;
                break;
            }
        }

        // 不考虑透射部分，走时已经超过当前最小走时，透射循环可结束
        if(sumt > travt) break;

        if(badrefrac) continue;

        // printf("down m=%d, refracted wave, sumt=%f, sumx=%f\n", m, sumt, sumx);

        // 统计走时 
        if(dist >= sumx){
            sumt += (dist - sumx)/v;
            if(sumt < travt){
                travt = sumt;
                // printf("refracted wave in layer %d, travt=%f\n", m, travt);
            }
        }

    } // END 寻找投射位置

    free(Vel);

    return travt;
}



/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt travt] %s\n\n", GRT_VERSION);printf(
"    A Supplementary Tool of GRT to Compute First Arrival Traveltime\n"
"    of P-wave and S-wave in Horizontally Layerd Halfspace Model. \n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt travt -M<model> -D<depsrc>/<deprcv> -R<r1>,<r2>[,...]\n"
"\n\n"
"Options:\n"
"----------------------------------------------------------------\n"
"    -M<model>    Filepath to 1D horizontally layered halfspace \n"
"                 model. The model file has 6 columns: \n"
"\n"
"         +-------+----------+----------+-------------+----+----+\n"
"         | H(km) | Vp(km/s) | Vs(km/s) | Rho(g/cm^3) | Qp | Qa |\n"
"         +-------+----------+----------+-------------+----+----+\n"
"\n"
"                 and the number of layers are unlimited.\n"
"\n"
"    -D<depsrc>/<deprcv>\n"
"                 <depsrc>: source depth (km).\n"
"                 <deprcv>: receiver depth (km).\n"
"\n"
"    -R<r1>,<r2>[,...]\n"
"                 Multiple epicentral distance (km), \n"
"                 seperated by comma.\n"
"\n"
"    -h           Display this help message.\n"
"\n\n"
"Examples:\n"
"----------------------------------------------------------------\n"
"    grt travt -Mmilrow -D2/0 -R10,20,30,40,50\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    char* command = Ctrl->name;
    int opt;
    while ((opt = getopt(argc, argv, ":M:D:R:h")) != -1) {
        switch (opt) {
            // 模型路径，其中每行分别为 
            //      厚度(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  Qp   Qs
            // 互相用空格隔开即可
            case 'M':
                Ctrl->M.active = true;
                Ctrl->M.s_modelpath = strdup(optarg);
                break;

            // 震源和场点深度， -Ddepsrc/deprcv
            case 'D':
                Ctrl->D.active = true;
                Ctrl->D.s_depsrc = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                Ctrl->D.s_deprcv = (char*)malloc(sizeof(char)*(strlen(optarg)+1));
                if(2 != sscanf(optarg, "%[^/]/%s", Ctrl->D.s_depsrc, Ctrl->D.s_deprcv)){
                    GRTBadOptionError(command, D, "");
                };
                if(1 != sscanf(Ctrl->D.s_depsrc, "%lf", &Ctrl->D.depsrc)){
                    GRTBadOptionError(command, D, "");
                }
                if(1 != sscanf(Ctrl->D.s_deprcv, "%lf", &Ctrl->D.deprcv)){
                    GRTBadOptionError(command, D, "");
                }
                if(Ctrl->D.depsrc < 0.0 || Ctrl->D.deprcv < 0.0){
                    GRTBadOptionError(command, D, "Negative value in -D is not supported.");
                }
                break;

            // 震中距数组，-Rr1,r2,r3,r4 ...
            case 'R':
                Ctrl->R.active = true;
                // 如果输入仅由数字、小数点和间隔符组成，则直接读取
                if(grt_string_composed_of(optarg, GRT_NUM_STR ".,")){
                    Ctrl->R.s_rs = grt_string_split(optarg, ",", &Ctrl->R.nr);
                } 
                // 否则从文件读取
                else {
                    FILE *fp = GRTCheckOpenFile(command, optarg, "r");
                    Ctrl->R.s_rs = grt_string_from_file(fp, &Ctrl->R.nr);
                    fclose(fp);
                }
                // 转为浮点数
                Ctrl->R.rs = (MYREAL*)realloc(Ctrl->R.rs, sizeof(MYREAL)*(Ctrl->R.nr));
                for(MYINT i=0; i<Ctrl->R.nr; ++i){
                    Ctrl->R.rs[i] = atof(Ctrl->R.s_rs[i]);
                    if(Ctrl->R.rs[i] < 0.0){
                        GRTBadOptionError(command, R, "Can't set negative epicentral distance(%f).", Ctrl->R.rs[i]);
                    }
                }
                break;

            GRT_Common_Options_in_Switch(command, (char)(optopt));
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(command, argc > 1);
    GRTCheckOptionActive(command, Ctrl, M);
    GRTCheckOptionActive(command, Ctrl, D);
    GRTCheckOptionActive(command, Ctrl, R);

}


/** 子模块主函数 */
int travt_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);

    getopt_from_command(Ctrl, argc, argv);

    // 读入模型文件
    if((Ctrl->M.mod1d = grt_read_mod1d_from_file(Ctrl->name, Ctrl->M.s_modelpath, Ctrl->D.depsrc, Ctrl->D.deprcv, true)) == NULL){
        exit(EXIT_FAILURE);
    }
    GRT_MODEL1D *mod1d = Ctrl->M.mod1d;

    printf("------------------------------------------------\n");
    printf(" Distance(km)     Tp(secs)         Ts(secs)     \n");
    double travtP=-1, travtS=-1;
    for(int i=0; i<Ctrl->R.nr; ++i){
        travtP = grt_compute_travt1d(
        mod1d->Thk, mod1d->Va, mod1d->n, mod1d->isrc, mod1d->ircv, Ctrl->R.rs[i]);
        travtS = grt_compute_travt1d(
        mod1d->Thk, mod1d->Vb, mod1d->n, mod1d->isrc, mod1d->ircv, Ctrl->R.rs[i]);
        
        printf(" %-15s  %-15.3f  %-15.3f\n", Ctrl->R.s_rs[i], travtP, travtS);
    }
    printf("------------------------------------------------\n");

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}