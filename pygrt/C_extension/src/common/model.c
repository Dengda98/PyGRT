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
    X(Thk, real_t)\
    X(Dep, real_t)\
    X(Va, real_t)\
    X(Vb, real_t)\
    X(Rho, real_t)\
    X(Qa, real_t)\
    X(Qb, real_t)\
    X(Qainv, real_t)\
    X(Qbinv, real_t)\
    X(mu, cplx_t)\
    X(lambda, cplx_t)\
    X(delta, cplx_t)\
    X(atna, cplx_t)\
    X(atnb, cplx_t)\
    X(xa, cplx_t)\
    X(xb, cplx_t)\
    X(caca, cplx_t)\
    X(cbcb, cplx_t)\


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
    for(size_t i=0; i<mod1d->n; ++i){
        if(i==mod1d->isrc){
            snprintf(indexstr, sizeof(indexstr), "%zu [src]", i+1);
        } else if(i==mod1d->ircv){
            snprintf(indexstr, sizeof(indexstr), "%zu [rcv]", i+1);
        } else {
            snprintf(indexstr, sizeof(indexstr), "%zu      ", i+1);
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


GRT_MODEL1D * grt_init_mod1d(size_t n){
    GRT_MODEL1D *mod1d = (GRT_MODEL1D *)calloc(1, sizeof(GRT_MODEL1D));
    mod1d->n = n;

    #define X(P, T)  mod1d->P = (T*)calloc(n, sizeof(T));
        GRT_FOR_EACH_MODEL_QUANTITY_ARRAY
    #undef X

    return mod1d;
}


GRT_MODEL1D * grt_copy_mod1d(const GRT_MODEL1D *mod1d1){
    GRT_MODEL1D *mod1d2 = (GRT_MODEL1D *)calloc(1, sizeof(GRT_MODEL1D));

    // 先直接赋值，实现浅拷贝
    *mod1d2 = *mod1d1;

    // 对指针部分再重新申请内存并赋值，实现深拷贝
    size_t n = mod1d1->n;
    #define X(P, T)  \
        mod1d2->P = (T*)calloc(n, sizeof(T));\
        memcpy(mod1d2->P, mod1d1->P, sizeof(T)*n);\

        GRT_FOR_EACH_MODEL_QUANTITY_ARRAY
    #undef X

    return mod1d2;
}


void grt_attenuate_mod1d(GRT_MODEL1D *mod1d, cplx_t omega){
    real_t Va0, Vb0;
    cplx_t atna, atnb;
    for(size_t i=0; i<mod1d->n; ++i){
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


void grt_mod1d_xa_xb(GRT_MODEL1D *mod1d, const real_t k)
{
    mod1d->k = k;
    // 不合理的频率值，只可能是在计算静态解，此时不需要xa, xb等物理量
    if(creal(mod1d->omega) < 0.0)  return;

    mod1d->c_phase = mod1d->omega/k;

    size_t isrc = mod1d->isrc;
    size_t ircv = mod1d->ircv;

    for(size_t i=0; i<mod1d->n; ++i){
        if( i == isrc || i == ircv ){
            mod1d->xa[i] = mod1d->xa[i-1];
            mod1d->caca[i] = mod1d->caca[i-1];
            mod1d->xb[i] = mod1d->xb[i-1];
            mod1d->cbcb[i] = mod1d->cbcb[i-1];
            continue;
        }

        real_t va, vb;
        va = mod1d->Va[i];
        vb = mod1d->Vb[i];
        cplx_t atna, atnb;
        atna = mod1d->atna[i];
        atnb = mod1d->atnb[i];

        cplx_t caca, cbcb;
        caca = mod1d->c_phase / (va*atna); 
        caca *= caca;
        mod1d->caca[i] = caca;
        mod1d->xa[i] = sqrt(1.0 - caca);
        
        cbcb = (vb > 0.0)? mod1d->c_phase / (vb*atnb) : 0.0;  // 考虑液体层
        cbcb *= cbcb;
        mod1d->cbcb[i] = cbcb;
        mod1d->xb[i] = sqrt(1.0 - cbcb);
    }
}


void grt_realloc_mod1d(GRT_MODEL1D *mod1d, size_t n){
    mod1d->n = n;

    #define X(P, T)  mod1d->P = (T*)realloc(mod1d->P, n*sizeof(T));
        GRT_FOR_EACH_MODEL_QUANTITY_ARRAY
    #undef X
}



GRT_MODEL1D * grt_read_mod1d_from_file(const char *modelpath, real_t depsrc, real_t deprcv, bool allowLiquid){
    GRTCheckFileExist(modelpath);
    
    FILE *fp = GRTCheckOpenFile(modelpath, "r");

    
    // 初始化
    GRT_MODEL1D *mod1d = grt_init_mod1d(1);

    const int ncols = 6; // 模型文件有6列，或除去qa qb有四列
    const int ncols_noQ = 4;
    size_t iline = 0;
    real_t h, va, vb, rho, qa, qb;
    real_t (*modarr)[ncols] = NULL;
    h = va = vb = rho = qa = qb = 0.0;
    size_t nlay = 0;
    mod1d->io_depth = false;

    size_t len;
    char *line = NULL;

    while(grt_getline(&line, &len, fp) != -1) {
        iline++;
        
        // 注释行
        if(grt_is_comment_or_empty(line))  continue;

        h = va = vb = rho = qa = qb = 0.0;
        int nscan = sscanf(line, "%lf %lf %lf %lf %lf %lf\n", &h, &va, &vb, &rho, &qa, &qb);
        if(ncols != nscan && ncols_noQ != nscan){
            GRTRaiseError("Model file read error in line %zu.\n", iline);
        };

        // 读取首行，如果首行首列为 0 ，则首列指示每层顶界面深度而非厚度
        if(nlay == 0 && h == 0.0){
            mod1d->io_depth = true;
        }

        if(va <= 0.0 || rho <= 0.0 || (ncols == nscan && (qa <= 0.0 || qb <= 0.0))){
            GRTRaiseError("In model file, line %zu, nonpositive value is not supported.\n", iline);
        }

        if(vb < 0.0){
            GRTRaiseError("In model file, line %zu, negative Vs is not supported.\n", iline);
        }

        if(!allowLiquid && vb == 0.0){
            GRTRaiseError("In model file, line %zu, Vs==0.0 is not supported.\n", iline);
        }

        modarr = (real_t(*)[ncols])realloc(modarr, sizeof(real_t)*ncols*(nlay+1));

        modarr[nlay][0] = h;
        modarr[nlay][1] = va;
        modarr[nlay][2] = vb;
        modarr[nlay][3] = rho;
        modarr[nlay][4] = qa;
        modarr[nlay][5] = qb;
        nlay++;

    }

    if(iline==0 || modarr==NULL){
        GRTRaiseError("Model file %s read error.\n", modelpath);
    }

    // 如果读取了深度，转为厚度
    if(mod1d->io_depth){
        for(size_t i=1; i<nlay; ++i){
            // 检查，若为负数，则表示输入的层顶深度非递增
            real_t tmp = modarr[i][0] - modarr[i-1][0];
            if(tmp < 0.0){
                GRTRaiseError("In model file, negative thickness found in layer %zu.\n", i);
            }
            modarr[i-1][0] = tmp;
        }
    }


    size_t isrc=0, ircv=0;
    size_t *pmin_idx, *pmax_idx, *pimg_idx;
    real_t depth = 0.0, depmin, depmax, depimg;
    bool ircvup = (depsrc >= deprcv);
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

    // 对最后一层的厚度做特殊处理
    modarr[nlay-1][0] = depmax + 1e30; // 保证够厚即可，用于下面定义虚拟层，实际计算不会用到最后一层厚度
    
    size_t nlay0 = nlay;
    nlay = 0;
    for(size_t i=0; i<nlay0; ++i){
        h = modarr[i][0];
        va = modarr[i][1];
        vb = modarr[i][2];
        rho = modarr[i][3];
        qa = modarr[i][4];
        qb = modarr[i][5];

        // 允许最后一层厚度为任意值
        if(h <= 0.0 && i < nlay0-1 ) {
            GRTRaiseError(
                "In line %zu, nonpositive thickness (except last layer) "
                "is not supported.\n", i+1);
        }

        // 划分震源层和接收层
        for(int k=0; k<2; ++k){
            // 这里的判断设计使得 min(isrc, ircv) >= 1，即不可能为 0 层
            // 这里不能改，否则会影响后续计算 R/T 矩阵的循环
            if(*pimg_idx == 0 && depth+h >= depimg && depsrc >= 0.0 && deprcv >= 0.0){
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
            "The receiver is located on the interface where there is liquid on one side. "
            "Due to the discontinuity of the tangential displacement on this interface, "
            "to reduce ambiguity, you should add a small offset to the receiver depth, "
            "thereby explicitly placing it within a specific layer. \n");
    }

    // 检查 --> 源点不能位于液-液、固-液界面
    if(isrc < nlay-1 && mod1d->Thk[isrc] == 0.0 && mod1d->Vb[isrc]*mod1d->Vb[isrc+1] == 0.0){
        GRTRaiseError(
            "The source is located on the interface where there is liquid on one side. "
            "Due to the discontinuity of the tangential displacement on this interface, "
            "to reduce ambiguity, you should add a small offset to the source depth, "
            "thereby explicitly placing it within a specific layer. \n");
    }

    // 将每层顶界面深度写入数组
    depth = 0.0;
    for(size_t iz=0; iz<mod1d->n; ++iz){
        mod1d->Dep[iz] = depth;
        depth += mod1d->Thk[iz];
    }

    fclose(fp);
    GRT_SAFE_FREE_PTR(modarr);
    GRT_SAFE_FREE_PTR(line);
    
    // 先指定负频率，仅填充弹性模量
    grt_attenuate_mod1d(mod1d, -1);

    return mod1d;
}


void grt_get_model_diglen_from_file(const char *modelpath, size_t diglen[6]){
    FILE *fp = GRTCheckOpenFile(modelpath, "r");
    size_t len;
    char *line = NULL;

    memset(diglen, 0, sizeof(size_t[6]));

    while(grt_getline(&line, &len, fp) != -1){
        char *token = strtok(line, " \n");
        for(int i=0; i<6; ++i){
            if(token == NULL) break;
            diglen[i] = GRT_MAX(diglen[i], strlen(token));
            token = strtok(NULL, " \n");
        }
    }

    GRT_SAFE_FREE_PTR(line);
    fclose(fp);
}


bool grt_check_vel_in_mod(const GRT_MODEL1D *mod1d, const real_t vel, const real_t tol){
    // 浮点数比较，检查是否存在该速度值
    for(size_t i=0; i<mod1d->n; ++i){
        if(fabs(vel - mod1d->Va[i])<tol || fabs(vel - mod1d->Vb[i])<tol)  return true;
    }
    return false;
}



void grt_get_mod1d_vmin_vmax(const GRT_MODEL1D *mod1d, real_t *vmin, real_t *vmax){
    *vmin = 9.0e30;
    *vmax = 0.0;
    const real_t *Va = mod1d->Va;
    const real_t *Vb = mod1d->Vb;
    for(size_t i=0; i<mod1d->n; ++i){
        if(Va[i] < *vmin) *vmin = Va[i];
        if(Va[i] > *vmax) *vmax = Va[i];
        if(Vb[i] < *vmin && Vb[i] > 0.0) *vmin = Vb[i];
        if(Vb[i] > *vmax && Vb[i] > 0.0) *vmax = Vb[i];
    }
}