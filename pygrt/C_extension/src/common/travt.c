/**
 * @file   travt.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-08
 * 
 *    计算一维均匀半无限层状介质的初至走时，思路借鉴了CPS330程序包的time96.f
 * 
 */

#include <stdbool.h>

#include "grt/common/travt.h"


real_t grt_compute_travt1d(
    const real_t *Thk, const real_t *Vel0, const size_t nlay, 
    const size_t isrc, const size_t ircv, const real_t dist)
{
    // 以防速度数组中存在零速度的情况，这里新建数组以去除0速度
    real_t *Vel = (real_t*)malloc(sizeof(real_t)*nlay);
    for(size_t i=0; i<nlay; ++i){
        Vel[i] = (Vel0[i] <= 0.0)? 1e-6 : Vel0[i];  // 给一个极慢值
    }

    // 根据互易原则，震源和场点可交换
    size_t imin, imax;
    imin = GRT_MIN(isrc, ircv);
    imax = GRT_MAX(isrc, ircv);

    // 如果震源与台站之间存在零速度层，则此时单一震相类型的走时不必再计算，
    // 直接返回默认值
    for(size_t i = imin; i <= imax; ++i){
        if(Vel0[i] == 0.0){
            GRT_SAFE_FREE_PTR(Vel);
            return -12345.00;
        }
    }
    

    // 震源和场点速度
    real_t vsrc = Vel[isrc];
    real_t vrcv = Vel[ircv];
    real_t vmax = (vsrc < vrcv)? vrcv : vsrc;
    real_t vmin = (vsrc < vrcv)? vsrc : vrcv;
    // 震源和场点深度
    real_t depsrc=0.0, deprcv=0.0;
    for(size_t i = 0; i < isrc; ++i) {depsrc += Thk[i];} 
    for(size_t i = 0; i < ircv; ++i) {deprcv += Thk[i];} 
    real_t depdif = fabs(depsrc - deprcv);


    // 初始化走时
    real_t travt = 9.0e30;

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
    if(imax - imin == 1){ // 位于同一物理层
        travt = sqrt(dist*dist + depdif*depdif) / vsrc;
        // printf("direct wave in same layer, travt=%f\n", travt);
    }
    else {
    // ------------------- 不同层直达波 ----------------------
    // -------------- 使用二分迭代法进行打靶 ------------------
        // 最大迭代次数
        const int nloop=50; 
        // 最小震中距差
        const real_t minX=1e-3;
        // 找到慢度上限，准确说是各层中最大慢度的最小值
        real_t pmax0=1.0/vmax;
        for(size_t i = imin; i <= imax; ++i){
            if(Thk[i] == 0.0) continue;
            if(pmax0 > 1.0/Vel[i])  pmax0 = 1.0/Vel[i];
        }
        // 初始化一些迭代变量
        real_t pmin=0.0, pmax=pmax0;
        real_t p;
        real_t s, c, v, h;
        real_t x = 0.0;
        real_t t = 0.0;
        real_t tint = 0.0;
        real_t dxdp = 0.0;
        for(int iter=0; iter<nloop; ++iter){
            x = t = tint = dxdp = 0.0;
            p = (pmin+pmax)/2.0;
            for(size_t i = imin; i < imax; ++i){
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
        real_t v, p, h, c;
        real_t sumt, sumx;
        bool badrefrac = false;
        // 找到透射位置
        for(size_t m = imin; m-- > 0;){ // 这样写是为了无符号整数的判断
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
            for(size_t i = imin; i < imax; ++i){
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
            for(size_t i = m+1; i < imin; ++i){
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

        }  // END 寻找透射位置

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
    for(size_t m = imax+1; m < nlay; ++m){
        real_t v, p, h, c;
        real_t sumt, sumx;
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
        for(size_t i = imin; i < imax; ++i){
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
        for(size_t i = imax; i < m; ++i){
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
