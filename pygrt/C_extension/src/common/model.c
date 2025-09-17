/**
 * @file   model.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * GRT_MODEL1D 结构体的相关操作函数
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "grt/common/model.h"
#include "grt/common/prtdbg.h"
#include "grt/common/attenuation.h"
#include "grt/common/util.h"

#include "grt/common/checkerror.h"


// 定义宏，方便写代码
#define GRT_FOR_EACH_MODEL_QUANTITY_ARRAY \
    X(Thk, MYREAL)\
    X(Dep, MYREAL)\
    X(Va, MYREAL)\
    X(Vb, MYREAL)\
    X(Rho, MYREAL)\
    X(Qa, MYREAL)\
    X(Qb, MYREAL)\
    X(Qainv, MYREAL)\
    X(Qbinv, MYREAL)\
    X(mu, MYCOMPLEX)\
    X(lambda, MYCOMPLEX)\
    X(delta, MYCOMPLEX)\
    X(atna, MYCOMPLEX)\
    X(atnb, MYCOMPLEX)\


void grt_print_mod1d(const GRT_MODEL1D *mod1d){
    // 模拟表格，打印速度
    // 每列字符宽度
    // [isrc/ircv] [h(km)] [Vp(km/s)] [Vs(km/s)] [Rho(g/cm^3)] [Qp] [Qs]
    const int ncols = 7;
    const int nlens[] = {13, 12, 13, 13, 16, 13, 13};
    int Nlen=0;
    for(int ic=0; ic<ncols; ++ic){
        Nlen += nlens[ic]; 
    }
    // 定义分割线
    char splitline[Nlen+2];
    {
        int n=0;
        for(int ic=0; ic<ncols; ++ic){
            splitline[n] = '+';
            for(int i=1; i<nlens[ic]; ++i){
                splitline[n + i] = '-';
            }
            n += nlens[ic];
        }
        splitline[Nlen] = '+';
        splitline[Nlen+1] = '\0';
    }
    printf("\n%s\n", splitline);

    // 打印题头
    printf("| %-*s ", nlens[0]-3, " ");
    printf("| %-*s ", nlens[1]-3, "H(km)");
    printf("| %-*s ", nlens[2]-3, "Vp(km/s)");
    printf("| %-*s ", nlens[3]-3, "Vs(km/s)");
    printf("| %-*s ", nlens[4]-3, "Rho(g/cm^3)");
    printf("| %-*s ", nlens[5]-3, "Qp");
    printf("| %-*s ", nlens[6]-3, "Qs");
    printf("|\n");
    printf("%s\n", splitline);


    char indexstr[nlens[0]-2+10];  // +10 以防止 -Wformat-truncation= 警告
    for(MYINT i=0; i<mod1d->n; ++i){
        if(i==mod1d->isrc){
            snprintf(indexstr, sizeof(indexstr), "%d [src]", i+1);
        } else if(i==mod1d->ircv){
            snprintf(indexstr, sizeof(indexstr), "%d [rcv]", i+1);
        } else {
            snprintf(indexstr, sizeof(indexstr), "%d      ", i+1);
        }

        printf("| %*s ", nlens[0]-3, indexstr);

        if(i < mod1d->n-1){
            printf("| %-*.2f ", nlens[1]-3, mod1d->Thk[i]);
        } else {
            printf("| %-*s ", nlens[1]-3, "Inf");
        }
        
        printf("| %-*.2f ", nlens[2]-3, mod1d->Va[i]);
        printf("| %-*.2f ", nlens[3]-3, mod1d->Vb[i]);
        printf("| %-*.2f ", nlens[4]-3, mod1d->Rho[i]);
        printf("| %-*.2e ", nlens[5]-3, mod1d->Qa[i]);
        printf("| %-*.2e ", nlens[6]-3, mod1d->Qb[i]);
        printf("|\n");
    }
    printf("%s\n", splitline);
    printf("\n");
}

void grt_free_mod1d(GRT_MODEL1D *mod1d){
    #define X(P, T)  GRT_SAFE_FREE_PTR(mod1d->P);
        GRT_FOR_EACH_MODEL_QUANTITY_ARRAY
    #undef X

    GRT_SAFE_FREE_PTR(mod1d);
}


GRT_MODEL1D * grt_init_mod1d(MYINT n){
    GRT_MODEL1D *mod1d = (GRT_MODEL1D *)calloc(1, sizeof(GRT_MODEL1D));
    mod1d->n = n;

    #define X(P, T)  mod1d->P = (T*)calloc(n, sizeof(T));
        GRT_FOR_EACH_MODEL_QUANTITY_ARRAY
    #undef X

    return mod1d;
}


GRT_MODEL1D * grt_copy_mod1d(const GRT_MODEL1D *mod1d1){
    GRT_MODEL1D *mod1d2 = grt_init_mod1d(mod1d1->n);
    MYINT n = mod1d1->n;
    mod1d2->n = mod1d1->n;
    mod1d2->depsrc = mod1d1->depsrc;
    mod1d2->deprcv = mod1d1->deprcv;
    mod1d2->isrc = mod1d1->isrc;
    mod1d2->ircv = mod1d1->ircv;
    mod1d2->ircvup = mod1d1->ircvup;
    mod1d2->io_depth = mod1d1->io_depth;

    #define X(P, T)  memcpy(mod1d2->P, mod1d1->P, sizeof(T)*n);
        GRT_FOR_EACH_MODEL_QUANTITY_ARRAY
    #undef X

    return mod1d2;
}


void grt_attenuate_mod1d(GRT_MODEL1D *mod1d, MYCOMPLEX omega){
    MYREAL Va0, Vb0;
    MYCOMPLEX atna, atnb;
    for(MYINT i=0; i<mod1d->n; ++i){
        Va0 = mod1d->Va[i];
        Vb0 = mod1d->Vb[i];

        // 圆频率实部为负数表明不考虑模型的 Q 值属性
        // 在读入模型后需要需要运行一次本函数以填充弹性模量，见 grt_read_mod1d_from_file 函数
        atna = (creal(omega) >= 0.0 && mod1d->Qainv[i] > 0.0)? grt_attenuation_law(mod1d->Qainv[i], omega) : 1.0;
        atnb = (creal(omega) >= 0.0 && mod1d->Qbinv[i] > 0.0)? grt_attenuation_law(mod1d->Qbinv[i], omega) : 1.0;

        mod1d->atna[i] = atna;
        mod1d->atnb[i] = atnb;
        
        mod1d->mu[i] = (Vb0*atnb)*(Vb0*atnb)*(mod1d->Rho[i]);
        mod1d->lambda[i] = (Va0*atnb)*(Va0*atnb)*(mod1d->Rho[i]) - 2*mod1d->mu[i];
        mod1d->delta[i] = (mod1d->lambda[i] + mod1d->mu[i]) / (mod1d->lambda[i] + 3.0*mod1d->mu[i]);
    }

#if Print_GRTCOEF == 1
    print_mod1d(mod1d);
#endif
}


void grt_get_mod1d_xa_xb(
    const GRT_MODEL1D *mod1d, const MYINT iy, const MYCOMPLEX c_phase, 
    MYCOMPLEX *pt_caca, MYCOMPLEX *pt_xa, MYCOMPLEX *pt_cbcb, MYCOMPLEX *pt_xb)
{
    MYREAL va, vb;
    va = mod1d->Va[iy];
    vb = mod1d->Vb[iy];
    MYCOMPLEX atna, atnb;
    atna = mod1d->atna[iy];
    atnb = mod1d->atnb[iy];

    MYCOMPLEX caca, cbcb;
    if(pt_caca!=NULL && pt_xa!=NULL){
        caca = c_phase / (va*atna); 
        caca *= caca;
        *pt_caca = caca;
        *pt_xa = sqrt(1.0 - caca);
    }
    
    if(pt_cbcb!=NULL && pt_xb!=NULL){
        cbcb = (vb > 0.0)? c_phase / (vb*atnb) : 0.0;  // 考虑液体层
        cbcb *= cbcb;
        *pt_cbcb = cbcb;
        *pt_xb = sqrt(1.0 - cbcb);
    }
}


void grt_realloc_mod1d(GRT_MODEL1D *mod1d, MYINT n){
    mod1d->n = n;

    #define X(P, T)  mod1d->P = (T*)realloc(mod1d->P, n*sizeof(T));
        GRT_FOR_EACH_MODEL_QUANTITY_ARRAY
    #undef X
}



GRT_MODEL1D * grt_read_mod1d_from_file(const char *command, const char *modelpath, double depsrc, double deprcv, bool allowLiquid){
    GRTCheckFileExist(command, modelpath);
    
    FILE *fp = GRTCheckOpenFile(command, modelpath, "r");

    MYINT isrc=-1, ircv=-1;
    MYINT *pmin_idx, *pmax_idx, *pimg_idx;
    double depth = 0.0, depmin, depmax, depimg;
    bool ircvup = (depsrc > deprcv);
    if(ircvup){
        pmin_idx = &ircv;
        pmax_idx = &isrc;
        depmin = deprcv;
        depmax = depsrc;
    } else {
        pmin_idx = &isrc;
        pmax_idx = &ircv;
        depmin = depsrc;
        depmax = deprcv;
    }
    depimg = depmin;
    pimg_idx = pmin_idx;

    // 初始化
    GRT_MODEL1D *mod1d = grt_init_mod1d(1);

    const int ncols = 6; // 模型文件有6列，或除去qa qb有四列
    const int ncols_noQ = 4;
    int iline = 0;
    double h, va, vb, rho, qa, qb;
    double (*modarr)[ncols] = NULL;
    h = va = vb = rho = qa = qb = 0.0;
    int nlay = 0;
    mod1d->io_depth = false;

    size_t len;
    char *line = NULL;

    while(grt_getline(&line, &len, fp) != -1) {
        iline++;
        
        // 注释行
        if(grt_is_comment_or_empty(line))  continue;

        h = va = vb = rho = qa = qb = 0.0;
        MYINT nscan = sscanf(line, "%lf %lf %lf %lf %lf %lf\n", &h, &va, &vb, &rho, &qa, &qb);
        if(ncols != nscan && ncols_noQ != nscan){
            GRTRaiseError("[%s] Model file read error in line %d.\n", command, iline);
        };

        // 读取首行，如果首行首列为 0 ，则首列指示每层顶界面深度而非厚度
        if(nlay == 0 && h == 0.0){
            mod1d->io_depth = true;
        }

        if(va <= 0.0 || rho <= 0.0 || (ncols == nscan && (qa <= 0.0 || qb <= 0.0))){
            GRTRaiseError("[%s] In model file, line %d, nonpositive value is not supported.\n", command, iline);
        }

        if(vb < 0.0){
            GRTRaiseError("[%s] In model file, line %d, negative Vs is not supported.\n", command, iline);
        }

        if(!allowLiquid && vb == 0.0){
            GRTRaiseError("[%s] In model file, line %d, Vs==0.0 is not supported.\n", command, iline);
        }

        modarr = (double(*)[ncols])realloc(modarr, sizeof(double)*ncols*(nlay+1));

        modarr[nlay][0] = h;
        modarr[nlay][1] = va;
        modarr[nlay][2] = vb;
        modarr[nlay][3] = rho;
        modarr[nlay][4] = qa;
        modarr[nlay][5] = qb;
        nlay++;

    }

    if(iline==0 || modarr==NULL){
        GRTRaiseError("[%s] Model file %s read error.\n", command, modelpath);
    }

    // 如果读取了深度，转为厚度
    if(mod1d->io_depth){
        for(int i=1; i<nlay; ++i){
            // 检查，若为负数，则表示输入的层顶深度非递增
            double tmp = modarr[i][0] - modarr[i-1][0];
            if(tmp < 0.0){
                GRTRaiseError("[%s] In model file, negative thickness found in layer %d.\n", command, i);
            }
            modarr[i-1][0] = tmp;
        }
    }

    // 对最后一层的厚度做特殊处理
    modarr[nlay-1][0] = depmax + 1e30; // 保证够厚即可，用于下面定义虚拟层，实际计算不会用到最后一层厚度
    
    int nlay0 = nlay;
    nlay = 0;
    for(int i=0; i<nlay0; ++i){
        h = modarr[i][0];
        va = modarr[i][1];
        vb = modarr[i][2];
        rho = modarr[i][3];
        qa = modarr[i][4];
        qb = modarr[i][5];

        // 允许最后一层厚度为任意值
        if(h <= 0.0 && i < nlay0-1 ) {
            GRTRaiseError("[%s] In line %d, nonpositive thickness (except last layer)"
                    " is not supported.\n", command, i+1);
        }

        // 划分震源层和接收层
        for(int k=0; k<2; ++k){
            // printf("%d, %d, %lf, %lf, %e ", i, k, depth+h, depimg, depth+h- depimg);
            if(*pimg_idx < 0 && depth+h >= depimg && depsrc >= 0.0 && deprcv >= 0.0){
                grt_realloc_mod1d(mod1d, nlay+1);
                mod1d->Thk[nlay] = depimg - depth;
                mod1d->Va[nlay] = va;
                mod1d->Vb[nlay] = vb;
                mod1d->Rho[nlay] = rho;
                mod1d->Qa[nlay] = qa;
                mod1d->Qb[nlay] = qb;
                mod1d->Qainv[nlay] = (qa > 0.0)? 1.0/qa : 0.0;
                mod1d->Qbinv[nlay] = (qb > 0.0)? 1.0/qb : 0.0;
                h = h - (depimg - depth);

                depth += depimg - depth;
                nlay++;

                depimg = depmax;
                *pimg_idx = nlay;
                pimg_idx = pmax_idx;
            }
        }
        

        grt_realloc_mod1d(mod1d, nlay+1);
        mod1d->Thk[nlay] = h;
        mod1d->Va[nlay] = va;
        mod1d->Vb[nlay] = vb;
        mod1d->Rho[nlay] = rho;
        mod1d->Qa[nlay] = qa;
        mod1d->Qb[nlay] = qb;
        mod1d->Qainv[nlay] = (qa > 0.0)? 1.0/qa : 0.0;
        mod1d->Qbinv[nlay] = (qb > 0.0)? 1.0/qb : 0.0;
        depth += h;
        nlay++;
    }

    mod1d->isrc = isrc;
    mod1d->ircv = ircv;
    mod1d->ircvup = ircvup;
    mod1d->n = nlay;
    mod1d->depsrc = depsrc;
    mod1d->deprcv = deprcv;

    // 检查，接收点不能位于液-液、固-液界面
    if(ircv < nlay-1 && mod1d->Thk[ircv] == 0.0 && mod1d->Vb[ircv]*mod1d->Vb[ircv+1] == 0.0){
        GRTRaiseError( 
            "[%s] The receiver is located on the interface where there is liquid on one side. "
            "Due to the discontinuity of the tangential displacement on this interface, "
            "to reduce ambiguity, you should add a small offset to the receiver depth, "
            "thereby explicitly placing it within a specific layer. \n", command);
    }

    // 检查 --> 源点不能位于液-液、固-液界面
    if(isrc < nlay-1 && mod1d->Thk[isrc] == 0.0 && mod1d->Vb[isrc]*mod1d->Vb[isrc+1] == 0.0){
        GRTRaiseError(
            "[%s] The source is located on the interface where there is liquid on one side. "
            "Due to the discontinuity of the tangential displacement on this interface, "
            "to reduce ambiguity, you should add a small offset to the source depth, "
            "thereby explicitly placing it within a specific layer. \n", command);
    }

    // 将每层顶界面深度写入数组
    depth = 0.0;
    for(int iz=0; iz<mod1d->n; ++iz){
        mod1d->Dep[iz] = depth;
        depth += mod1d->Thk[iz];
    }

    fclose(fp);
    GRT_SAFE_FREE_PTR(modarr);
    GRT_SAFE_FREE_PTR(line);
    
    // 填充弹性模量
    grt_attenuate_mod1d(mod1d, -1);

    return mod1d;
}


void grt_get_model_diglen_from_file(const char *command, const char *modelpath, MYINT diglen[6]){
    FILE *fp = GRTCheckOpenFile(command, modelpath, "r");
    size_t len;
    char *line = NULL;

    memset(diglen, 0, sizeof(MYINT)*6);

    while(grt_getline(&line, &len, fp) != -1){
        char *token = strtok(line, " \n");
        for(MYINT i=0; i<6; ++i){
            if(token == NULL) break;
            diglen[i] = GRT_MAX(diglen[i], (MYINT)strlen(token));
            token = strtok(NULL, " \n");
        }
    }

    GRT_SAFE_FREE_PTR(line);
    fclose(fp);
}


bool grt_check_vel_in_mod(const GRT_MODEL1D *mod1d, const MYREAL vel, const MYREAL tol){
    // 浮点数比较，检查是否存在该速度值
    for(MYINT i=0; i<mod1d->n; ++i){
        if(fabs(vel - mod1d->Va[i])<tol || fabs(vel - mod1d->Vb[i])<tol)  return true;
    }
    return false;
}



void grt_get_mod1d_vmin_vmax(const GRT_MODEL1D *mod1d, MYREAL *vmin, MYREAL *vmax){
    *vmin = 9.0e30;
    *vmax = 0.0;
    const MYREAL *Va = mod1d->Va;
    const MYREAL *Vb = mod1d->Vb;
    for(MYINT i=0; i<mod1d->n; ++i){
        if(Va[i] < *vmin) *vmin = Va[i];
        if(Va[i] > *vmax) *vmax = Va[i];
        if(Vb[i] < *vmin && Vb[i] > 0.0) *vmin = Vb[i];
        if(Vb[i] > *vmax && Vb[i] > 0.0) *vmax = Vb[i];
    }
}