/**
 * @file   model.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * MODEL1D结构体的相关操作函数
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <complex.h>

#include "grt/common/model.h"
#include "grt/common/prtdbg.h"
#include "grt/common/attenuation.h"
#include "grt/common/colorstr.h"

#include "grt/common/checkerror.h"

void grt_print_mod1d(const GRT_MODEL1D *mod1d){
    LAYER *lay;
    for(MYINT u=0; u<50; ++u){printf("---"); } printf("\n");
    for(MYINT i=0; i<mod1d->n; ++i){
        lay = mod1d->lays+i;
        printf("     Dep=%6.2f, Va=%6.2f, Vb=%6.2f, thk=%6.2f, Rho=%6.2f, 1/Qa=%6.2e, 1/Qb=%6.2e\n",
                     lay->dep, lay->Va, lay->Vb, lay->thk, lay->Rho, lay->Qainv, lay->Qbinv);
        printf("     mu=(%e %+e I)\n", creal(lay->mu), cimag(lay->mu));
        printf("     lambda=(%e %+e I)\n", creal(lay->lambda), cimag(lay->lambda));
        printf("     delta=(%e %+e I)\n", creal(lay->delta), cimag(lay->delta));
        printf("     ka^2=%e%+eJ\n", creal(lay->kaka), cimag(lay->kaka));
        printf("     kb^2=%e%+eJ\n", creal(lay->kbkb), cimag(lay->kbkb));
        for(MYINT u=0; u<50; ++u){printf("---"); } printf("\n");
    }
}

void grt_print_pymod(const GRT_PYMODEL1D *pymod){
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
    for(MYINT i=0; i<pymod->n; ++i){
        if(i==pymod->isrc){
            snprintf(indexstr, sizeof(indexstr), "%d [src]", i+1);
        } else if(i==pymod->ircv){
            snprintf(indexstr, sizeof(indexstr), "%d [rcv]", i+1);
        } else {
            snprintf(indexstr, sizeof(indexstr), "%d      ", i+1);
        }

        printf("| %*s ", nlens[0]-3, indexstr);

        if(i < pymod->n-1){
            printf("| %-*.2f ", nlens[1]-3, pymod->Thk[i]);
        } else {
            printf("| %-*s ", nlens[1]-3, "Inf");
        }
        
        printf("| %-*.2f ", nlens[2]-3, pymod->Va[i]);
        printf("| %-*.2f ", nlens[3]-3, pymod->Vb[i]);
        printf("| %-*.2f ", nlens[4]-3, pymod->Rho[i]);
        printf("| %-*.2e ", nlens[5]-3, pymod->Qa[i]);
        printf("| %-*.2e ", nlens[6]-3, pymod->Qb[i]);
        printf("|\n");
    }
    printf("%s\n", splitline);
    printf("\n");
}

void grt_free_pymod(GRT_PYMODEL1D *pymod){
    GRT_SAFE_FREE_PTR(pymod->Thk);
    GRT_SAFE_FREE_PTR(pymod->Va);
    GRT_SAFE_FREE_PTR(pymod->Vb);
    GRT_SAFE_FREE_PTR(pymod->Rho);
    GRT_SAFE_FREE_PTR(pymod->Qa);
    GRT_SAFE_FREE_PTR(pymod->Qb);
    GRT_SAFE_FREE_PTR(pymod);
}


GRT_MODEL1D * grt_init_mod1d(MYINT n){
    GRT_MODEL1D *mod1d = (GRT_MODEL1D *)malloc(sizeof(GRT_MODEL1D));
    mod1d->n = n;
    mod1d->lays = (LAYER *)malloc(sizeof(LAYER) * n);
    return mod1d;
}



void grt_get_mod1d(const GRT_PYMODEL1D *pymod1d, GRT_MODEL1D *mod1d){
    MYINT n = pymod1d->n;
    mod1d->n = n;
    mod1d->isrc = pymod1d->isrc;
    mod1d->ircv = pymod1d->ircv;
    mod1d->ircvup = pymod1d->ircvup;

    if(mod1d->ircvup){
        mod1d->imin = mod1d->ircv; 
        mod1d->imax = mod1d->isrc; 
    } else {
        mod1d->imin = mod1d->isrc; 
        mod1d->imax = mod1d->ircv; 
    }

    // MYREAL Rho0;
    // MYREAL Vb0;
    LAYER *lay;
    MYREAL dep=0.0;
    for(MYINT i=0; i<n; ++i){
        lay = mod1d->lays + i;
        lay->thk = pymod1d->Thk[i];
        lay->dep = dep;
        lay->Va  = pymod1d->Va[i];
        lay->Vb  = pymod1d->Vb[i];
        lay->Rho = pymod1d->Rho[i];
        lay->Qainv  = (pymod1d->Qa[i] > 0.0)? RONE/pymod1d->Qa[i] : 0.0;
        lay->Qbinv  = (pymod1d->Qb[i] > 0.0)? RONE/pymod1d->Qb[i] : 0.0;

        lay->mu = (lay->Vb)*(lay->Vb)*(lay->Rho);
        lay->lambda = (lay->Va)*(lay->Va)*(lay->Rho) - RTWO*lay->mu;
        lay->delta = (lay->lambda + lay->mu) / (lay->lambda + RTHREE*lay->mu);

        dep += pymod1d->Thk[i];
    }

}


void grt_copy_mod1d(const GRT_MODEL1D *mod1d1, GRT_MODEL1D *mod1d2){
    MYINT n = mod1d1->n;
    mod1d2->n = mod1d1->n;
    mod1d2->isrc = mod1d1->isrc;
    mod1d2->ircv = mod1d1->ircv;
    mod1d2->ircvup = mod1d1->ircvup;

    mod1d2->imin = mod1d1->imin;
    mod1d2->imax = mod1d1->imax;

    LAYER *lay1, *lay2;
    for(MYINT i=0; i<n; ++i){
        lay1 = mod1d1->lays + i;
        lay2 = mod1d2->lays + i;

        lay2->thk = lay1->thk;
        lay2->dep = lay1->dep;
        lay2->Va = lay1->Va;
        lay2->Vb = lay1->Vb;
        lay2->Rho = lay1->Rho;
        lay2->Qainv = lay1->Qainv;
        lay2->Qbinv = lay1->Qbinv;

        lay2->mu = lay1->mu;
        lay2->kaka = lay1->kaka;
        lay2->kbkb = lay1->kbkb;

        lay2->lambda = lay1->lambda;
        lay2->delta = lay1->delta;
    }
    
}


void grt_free_mod1d(GRT_MODEL1D *mod1d){
    GRT_SAFE_FREE_PTR(mod1d->lays);
    GRT_SAFE_FREE_PTR(mod1d);
}



void grt_update_mod1d_omega(GRT_MODEL1D *mod1d, MYCOMPLEX omega){
    MYREAL Va0, Vb0;
    MYCOMPLEX ka0, kb0;
    MYCOMPLEX atna, atnb;
    LAYER *lay;
    for(MYINT i=0; i<mod1d->n; ++i){
        lay = mod1d->lays + i;
        Va0 = lay->Va;
        Vb0 = lay->Vb;

        atna = (lay->Qainv > 0.0)? grt_attenuation_law(lay->Qainv, omega) : 1.0;
        atnb = (lay->Qbinv > 0.0)? grt_attenuation_law(lay->Qbinv, omega) : 1.0;
        
        ka0 = omega/(Va0*atna);
        kb0 = (Vb0>RZERO)? omega/(Vb0*atnb) : CZERO;
        lay->kaka = ka0*ka0;
        lay->kbkb = kb0*kb0;
        
        lay->mu = (Vb0*atnb)*(Vb0*atnb)*(lay->Rho);
        lay->lambda = (Va0*atnb)*(Va0*atnb)*(lay->Rho) - 2*lay->mu;
        lay->delta = (lay->lambda + lay->mu) / (lay->lambda + 3*lay->mu);
    }

#if Print_GRTCOEF == 1
    print_mod1d(mod1d);
#endif
}


GRT_PYMODEL1D * grt_init_pymod(MYINT n){
    GRT_PYMODEL1D *pymod = (GRT_PYMODEL1D *)calloc(1, sizeof(GRT_PYMODEL1D));
    pymod->n = n;
    
    pymod->Thk = (MYREAL*)calloc(n, sizeof(MYREAL));
    pymod->Va = (MYREAL*)calloc(n, sizeof(MYREAL));
    pymod->Vb = (MYREAL*)calloc(n, sizeof(MYREAL));
    pymod->Rho = (MYREAL*)calloc(n, sizeof(MYREAL));
    pymod->Qa = (MYREAL*)calloc(n, sizeof(MYREAL));
    pymod->Qb = (MYREAL*)calloc(n, sizeof(MYREAL));

    return pymod;
}

void grt_realloc_pymod(GRT_PYMODEL1D *pymod, MYINT n){
    pymod->n = n;

    pymod->Thk = (MYREAL*)realloc(pymod->Thk, n*sizeof(MYREAL));
    pymod->Va = (MYREAL*)realloc(pymod->Va, n*sizeof(MYREAL));
    pymod->Vb = (MYREAL*)realloc(pymod->Vb, n*sizeof(MYREAL));
    pymod->Rho = (MYREAL*)realloc(pymod->Rho, n*sizeof(MYREAL));
    pymod->Qa = (MYREAL*)realloc(pymod->Qa, n*sizeof(MYREAL));
    pymod->Qb = (MYREAL*)realloc(pymod->Qb, n*sizeof(MYREAL));
}


GRT_PYMODEL1D * grt_read_pymod_from_file(const char *command, const char *modelpath, double depsrc, double deprcv, bool allowLiquid){
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
    GRT_PYMODEL1D *pymod = grt_init_pymod(1);

    const int ncols = 6; // 模型文件有6列，或除去qa qb有四列
    const int ncols_noQ = 4;
    char line[1024];
    int iline = 0;
    double h, va, vb, rho, qa, qb;
    double (*modarr)[ncols] = NULL;
    h = va = vb = rho = qa = qb = 0.0;
    int nlay = 0;
    pymod->io_depth = false;

    while(fgets(line, sizeof(line), fp)) {
        iline++;
        
        // 注释行
        if(line[0]=='#')  continue;

        h = va = vb = rho = qa = qb = 0.0;
        MYINT nscan = sscanf(line, "%lf %lf %lf %lf %lf %lf\n", &h, &va, &vb, &rho, &qa, &qb);
        if(ncols != nscan && ncols_noQ != nscan){
            fprintf(stderr, "[%s] " BOLD_RED "Model file read error in line %d.\n" DEFAULT_RESTORE, command, iline);
            return NULL;
        };

        // 读取首行，如果首行首列为 0 ，则首列指示每层顶界面深度而非厚度
        if(nlay == 0 && h == 0.0){
            pymod->io_depth = true;
        }

        if(va <= 0.0 || rho <= 0.0 || (ncols == nscan && (qa <= 0.0 || qb <= 0.0))){
            fprintf(stderr, "[%s] " BOLD_RED "In model file, line %d, nonpositive value is not supported.\n" DEFAULT_RESTORE, command, iline);
            return NULL;
        }

        if(vb < 0.0){
            fprintf(stderr, "[%s] " BOLD_RED "In model file, line %d, negative Vs is not supported.\n" DEFAULT_RESTORE, command, iline);
            return NULL;
        }

        if(!allowLiquid && vb == 0.0){
            fprintf(stderr, "[%s] " BOLD_RED "In model file, line %d, Vs==0.0 is not supported.\n" DEFAULT_RESTORE, command, iline);
            return NULL;
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
        fprintf(stderr, "[%s] " BOLD_RED "Model file read error.\n" DEFAULT_RESTORE, command);
        return NULL;
    }

    // 如果读取了深度，转为厚度
    if(pymod->io_depth){
        for(int i=1; i<nlay; ++i){
            modarr[i-1][0] = modarr[i][0] - modarr[i-1][0];
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
            fprintf(stderr, "[%s] " BOLD_RED "In line %d, nonpositive thickness (except last layer)"
                    " is not supported.\n" DEFAULT_RESTORE, command, i+1);
            return NULL;
        }

        // 划分震源层和接收层
        for(int k=0; k<2; ++k){
            // printf("%d, %d, %lf, %lf, %e ", i, k, depth+h, depimg, depth+h- depimg);
            if(*pimg_idx < 0 && depth+h >= depimg && depsrc >= 0.0 && deprcv >= 0.0){
                grt_realloc_pymod(pymod, nlay+1);
                pymod->Thk[nlay] = depimg - depth;
                pymod->Va[nlay] = va;
                pymod->Vb[nlay] = vb;
                pymod->Rho[nlay] = rho;
                pymod->Qa[nlay] = qa;
                pymod->Qb[nlay] = qb;
                h = h - (depimg - depth);

                depth += depimg - depth;
                nlay++;

                depimg = depmax;
                *pimg_idx = nlay;
                pimg_idx = pmax_idx;
            }
        }
        

        grt_realloc_pymod(pymod, nlay+1);
        pymod->Thk[nlay] = h;
        pymod->Va[nlay] = va;
        pymod->Vb[nlay] = vb;
        pymod->Rho[nlay] = rho;
        pymod->Qa[nlay] = qa;
        pymod->Qb[nlay] = qb;
        depth += h;
        nlay++;
    }

    pymod->isrc = isrc;
    pymod->ircv = ircv;
    pymod->ircvup = ircvup;
    pymod->n = nlay;
    pymod->depsrc = depsrc;
    pymod->deprcv = deprcv;

    // 检查，接收点不能位于液-液、固-液界面
    if(ircv < nlay-1 && pymod->Thk[ircv] == 0.0 && pymod->Vb[ircv]*pymod->Vb[ircv+1] == 0.0){
        fprintf(stderr, 
            "[%s] " BOLD_RED "The receiver is located on the interface where there is liquid on one side. "
            "Due to the discontinuity of the tangential displacement on this interface, "
            "to reduce ambiguity, you should add a small offset to the receiver depth, "
            "thereby explicitly placing it within a specific layer. \n"
            DEFAULT_RESTORE, command);
        return NULL;
    }

    // 检查 --> 源点不能位于液-液、固-液界面
    if(isrc < nlay-1 && pymod->Thk[isrc] == 0.0 && pymod->Vb[isrc]*pymod->Vb[isrc+1] == 0.0){
        fprintf(stderr, 
            "[%s] " BOLD_RED "The source is located on the interface where there is liquid on one side. "
            "Due to the discontinuity of the tangential displacement on this interface, "
            "to reduce ambiguity, you should add a small offset to the source depth, "
            "thereby explicitly placing it within a specific layer. \n"
            DEFAULT_RESTORE, command);
        return NULL;
    }

    fclose(fp);
    GRT_SAFE_FREE_PTR(modarr);
    
    return pymod;
}


void grt_get_pymod_vmin_vmax(const GRT_PYMODEL1D *pymod, double *vmin, double *vmax){
    *vmin = __DBL_MAX__;
    *vmax = RZERO;
    const MYREAL *Va = pymod->Va;
    const MYREAL *Vb = pymod->Vb;
    for(MYINT i=0; i<pymod->n; ++i){
        if(Va[i] < *vmin) *vmin = Va[i];
        if(Va[i] > *vmax) *vmax = Va[i];
        if(Vb[i] < *vmin && Vb[i] > RZERO) *vmin = Vb[i];
        if(Vb[i] > *vmax && Vb[i] > RZERO) *vmax = Vb[i];
    }
}