/**
 * @file   stgrnlib.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态格林函数库 STGRNLIB：内存管理与四维 nc 读写
 *
 */

#include <string.h>

#include "grt/static/stgrnlib.h"
#include "grt/common/mynetcdf.h"
#include "grt/common/checkerror.h"

/** 按已填好的维度申请 u/uiz/uir */
static void allocate_u(STGRNLIB *lib)
{
    size_t nr = lib->nnorth * lib->neast;
    lib->u = (realChnlGrid ***)calloc(lib->ndepsrc, sizeof(*lib->u));
    lib->uiz = lib->calc_upar ? (realChnlGrid ***)calloc(lib->ndepsrc, sizeof(*lib->uiz)) : NULL;
    lib->uir = lib->calc_upar ? (realChnlGrid ***)calloc(lib->ndepsrc, sizeof(*lib->uir)) : NULL;

    for(size_t is = 0; is < lib->ndepsrc; ++is){
        lib->u[is] = (realChnlGrid **)calloc(lib->ndeprcv, sizeof(*lib->u[is]));
        if(lib->calc_upar){
            lib->uiz[is] = (realChnlGrid **)calloc(lib->ndeprcv, sizeof(*lib->uiz[is]));
            lib->uir[is] = (realChnlGrid **)calloc(lib->ndeprcv, sizeof(*lib->uir[is]));
        }
        for(size_t ir = 0; ir < lib->ndeprcv; ++ir){
            lib->u[is][ir] = (realChnlGrid *)calloc(nr, sizeof(realChnlGrid));
            if(lib->calc_upar){
                lib->uiz[is][ir] = (realChnlGrid *)calloc(nr, sizeof(realChnlGrid));
                lib->uir[is][ir] = (realChnlGrid *)calloc(nr, sizeof(realChnlGrid));
            }
        }
    }
}

/** 仅释放 u/uiz/uir */
static void free_u(STGRNLIB *lib)
{
    if(lib == NULL) return;
    if(lib->u != NULL){
        for(size_t is = 0; is < lib->ndepsrc; ++is){
            if(lib->u[is] != NULL){
                for(size_t ir = 0; ir < lib->ndeprcv; ++ir){
                    GRT_SAFE_FREE_PTR(lib->u[is][ir]);
                    if(lib->uiz && lib->uiz[is]) GRT_SAFE_FREE_PTR(lib->uiz[is][ir]);
                    if(lib->uir && lib->uir[is]) GRT_SAFE_FREE_PTR(lib->uir[is][ir]);
                }
            }
            GRT_SAFE_FREE_PTR(lib->u[is]);
            if(lib->uiz) GRT_SAFE_FREE_PTR(lib->uiz[is]);
            if(lib->uir) GRT_SAFE_FREE_PTR(lib->uir[is]);
        }
    }
    GRT_SAFE_FREE_PTR(lib->u);
    GRT_SAFE_FREE_PTR(lib->uiz);
    GRT_SAFE_FREE_PTR(lib->uir);
}

STGRNLIB *grt_stgrnlib_alloc(
    size_t ndepsrc, const real_t *depsrcs,
    size_t ndeprcv, const real_t *deprcvs,
    size_t nnorth,  const real_t *norths,
    size_t neast,   const real_t *easts,
    bool calc_upar)
{
    if(ndepsrc == 0 || ndeprcv == 0 || nnorth == 0 || neast == 0){
        GRTRaiseError("dimensions must be positive.");
    }
    if(depsrcs == NULL || deprcvs == NULL || norths == NULL || easts == NULL){
        GRTRaiseError("coordinate arrays are NULL.");
    }

    STGRNLIB *lib = (STGRNLIB *)calloc(1, sizeof(*lib));
    lib->ndepsrc = ndepsrc;
    lib->ndeprcv = ndeprcv;
    lib->nnorth = nnorth;
    lib->neast = neast;
    lib->calc_upar = calc_upar;

    lib->depsrcs = (real_t *)malloc(ndepsrc * sizeof(real_t));
    lib->deprcvs = (real_t *)malloc(ndeprcv * sizeof(real_t));
    lib->norths = (real_t *)malloc(nnorth * sizeof(real_t));
    lib->easts = (real_t *)malloc(neast * sizeof(real_t));
    memcpy(lib->depsrcs, depsrcs, ndepsrc * sizeof(real_t));
    memcpy(lib->deprcvs, deprcvs, ndeprcv * sizeof(real_t));
    memcpy(lib->norths, norths, nnorth * sizeof(real_t));
    memcpy(lib->easts, easts, neast * sizeof(real_t));

    // 介质参数由调用方随后填入
    lib->src_va = (real_t *)calloc(ndepsrc, sizeof(real_t));
    lib->src_vb = (real_t *)calloc(ndepsrc, sizeof(real_t));
    lib->src_rho = (real_t *)calloc(ndepsrc, sizeof(real_t));
    lib->rcv_va = (real_t *)calloc(ndeprcv, sizeof(real_t));
    lib->rcv_vb = (real_t *)calloc(ndeprcv, sizeof(real_t));
    lib->rcv_rho = (real_t *)calloc(ndeprcv, sizeof(real_t));

    allocate_u(lib);
    return lib;
}

void grt_stgrnlib_free(STGRNLIB *lib)
{
    if(lib == NULL) return;
    free_u(lib);
    GRT_SAFE_FREE_PTR(lib->depsrcs);
    GRT_SAFE_FREE_PTR(lib->deprcvs);
    GRT_SAFE_FREE_PTR(lib->norths);
    GRT_SAFE_FREE_PTR(lib->easts);
    GRT_SAFE_FREE_PTR(lib->src_va);
    GRT_SAFE_FREE_PTR(lib->src_vb);
    GRT_SAFE_FREE_PTR(lib->src_rho);
    GRT_SAFE_FREE_PTR(lib->rcv_va);
    GRT_SAFE_FREE_PTR(lib->rcv_vb);
    GRT_SAFE_FREE_PTR(lib->rcv_rho);
    free(lib);
}

/** 读入一层 (is, ir) 的通道；ncid 已打开，变量为 4D */
static void read_nc_channels_slice(STGRNLIB *lib, size_t is, size_t ir, int ncid)
{
    size_t nr = lib->nnorth * lib->neast;
    real_t *buf = (real_t *)calloc(nr, sizeof(real_t));
    size_t start[4] = {is, ir, 0, 0};
    size_t count[4] = {1, 1, lib->nnorth, lib->neast};

    GRT_LOOP_ChnlGrid(im, c){
        int modr = GRT_SRC_M_ORDERS[im];
        if(modr == 0 && GRT_ZRT_CODES[c] == 'T') continue;

        char *s_title = NULL;
        int varid;

        GRT_SAFE_ASPRINTF(&s_title, "%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
        NC_CHECK(nc_inq_varid(ncid, s_title, &varid));
        NC_CHECK(NC_FUNC_REAL(nc_get_vara)(ncid, varid, start, count, buf));
        for(size_t ipt = 0; ipt < nr; ++ipt){
            lib->u[is][ir][ipt][im][c] = buf[ipt];
        }
        GRT_SAFE_FREE_PTR(s_title);

        if(lib->calc_upar){
            GRT_SAFE_ASPRINTF(&s_title, "z%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
            NC_CHECK(nc_inq_varid(ncid, s_title, &varid));
            NC_CHECK(NC_FUNC_REAL(nc_get_vara)(ncid, varid, start, count, buf));
            for(size_t ipt = 0; ipt < nr; ++ipt){
                lib->uiz[is][ir][ipt][im][c] = buf[ipt];
            }
            GRT_SAFE_FREE_PTR(s_title);

            GRT_SAFE_ASPRINTF(&s_title, "r%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
            NC_CHECK(nc_inq_varid(ncid, s_title, &varid));
            NC_CHECK(NC_FUNC_REAL(nc_get_vara)(ncid, varid, start, count, buf));
            for(size_t ipt = 0; ipt < nr; ++ipt){
                lib->uir[is][ir][ipt][im][c] = buf[ipt];
            }
            GRT_SAFE_FREE_PTR(s_title);
        }
    }

    GRT_SAFE_FREE_PTR(buf);
}

STGRNLIB *grt_stgrnlib_load_nc(const char *path)
{
    GRTCheckFileExist(path);

    int ncid;
    NC_CHECK(nc_open(path, NC_NOWRITE, &ncid));

    int depsrc_dimid, deprcv_dimid, north_dimid, east_dimid;
    size_t ndepsrc, ndeprcv, nnorth, neast;
    NC_CHECK(nc_inq_dimid(ncid, "depsrc", &depsrc_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, depsrc_dimid, &ndepsrc));
    NC_CHECK(nc_inq_dimid(ncid, "deprcv", &deprcv_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, deprcv_dimid, &ndeprcv));
    NC_CHECK(nc_inq_dimid(ncid, "north", &north_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, north_dimid, &nnorth));
    NC_CHECK(nc_inq_dimid(ncid, "east", &east_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, east_dimid, &neast));

    if(ndepsrc == 0 || ndeprcv == 0 || nnorth == 0 || neast == 0){
        NC_CHECK(nc_close(ncid));
        GRTRaiseError("Invalid STGRNLIB nc \"%s\": empty dimension.", path);
    }

    int int_calc_upar = 0;
    NC_CHECK(nc_get_att_int(ncid, NC_GLOBAL, "calc_upar", &int_calc_upar));

    // 坐标轴读入后经 grt_stgrnlib_alloc 建壳
    real_t *depsrcs = (real_t *)calloc(ndepsrc, sizeof(real_t));
    real_t *deprcvs = (real_t *)calloc(ndeprcv, sizeof(real_t));
    real_t *norths = (real_t *)calloc(nnorth, sizeof(real_t));
    real_t *easts = (real_t *)calloc(neast, sizeof(real_t));

    int varid;
    NC_CHECK(nc_inq_varid(ncid, "depsrc", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, depsrcs));
    NC_CHECK(nc_inq_varid(ncid, "deprcv", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, deprcvs));
    NC_CHECK(nc_inq_varid(ncid, "north", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, norths));
    NC_CHECK(nc_inq_varid(ncid, "east", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, easts));

    STGRNLIB *lib = grt_stgrnlib_alloc(
        ndepsrc, depsrcs, ndeprcv, deprcvs, nnorth, norths, neast, easts, (int_calc_upar != 0));
    GRT_SAFE_FREE_PTR(depsrcs);
    GRT_SAFE_FREE_PTR(deprcvs);
    GRT_SAFE_FREE_PTR(norths);
    GRT_SAFE_FREE_PTR(easts);

    NC_CHECK(nc_inq_varid(ncid, "src_va", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, lib->src_va));
    NC_CHECK(nc_inq_varid(ncid, "src_vb", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, lib->src_vb));
    NC_CHECK(nc_inq_varid(ncid, "src_rho", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, lib->src_rho));
    NC_CHECK(nc_inq_varid(ncid, "rcv_va", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, lib->rcv_va));
    NC_CHECK(nc_inq_varid(ncid, "rcv_vb", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, lib->rcv_vb));
    NC_CHECK(nc_inq_varid(ncid, "rcv_rho", &varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, lib->rcv_rho));

    for(size_t is = 0; is < ndepsrc; ++is){
        for(size_t ir = 0; ir < ndeprcv; ++ir){
            read_nc_channels_slice(lib, is, ir, ncid);
        }
    }

    NC_CHECK(nc_close(ncid));
    return lib;
}

void grt_stgrnlib_save_nc(const STGRNLIB *lib, const char *path)
{
    if(lib == NULL || path == NULL){
        GRTRaiseError("lib/path is NULL.");
    }
    if(lib->ndepsrc == 0 || lib->ndeprcv == 0 || lib->nnorth == 0 || lib->neast == 0){
        GRTRaiseError("empty library.");
    }

    int ncid;
    int depsrc_dimid, deprcv_dimid, north_dimid, east_dimid;
    int dimids[4];
    int depsrc_varid, deprcv_varid, north_varid, east_varid;
    int src_va_varid, src_vb_varid, src_rho_varid;
    int rcv_va_varid, rcv_vb_varid, rcv_rho_varid;
    intChnlGrid u_varids;
    intChnlGrid uiz_varids;
    intChnlGrid uir_varids;

    NC_CHECK(nc_create(path, NC_CLOBBER, &ncid));

    {
        int tmp = lib->calc_upar ? 1 : 0;
        NC_CHECK(nc_put_att_int(ncid, NC_GLOBAL, "calc_upar", NC_INT, 1, &tmp));
    }

    NC_CHECK(nc_def_dim(ncid, "depsrc", lib->ndepsrc, &depsrc_dimid));
    NC_CHECK(nc_def_dim(ncid, "deprcv", lib->ndeprcv, &deprcv_dimid));
    NC_CHECK(nc_def_dim(ncid, "north", lib->nnorth, &north_dimid));
    NC_CHECK(nc_def_dim(ncid, "east", lib->neast, &east_dimid));
    dimids[0] = depsrc_dimid;
    dimids[1] = deprcv_dimid;
    dimids[2] = north_dimid;
    dimids[3] = east_dimid;

    NC_CHECK(nc_def_var(ncid, "depsrc", NC_REAL, 1, &depsrc_dimid, &depsrc_varid));
    NC_CHECK(nc_def_var(ncid, "deprcv", NC_REAL, 1, &deprcv_dimid, &deprcv_varid));
    NC_CHECK(nc_def_var(ncid, "north", NC_REAL, 1, &north_dimid, &north_varid));
    NC_CHECK(nc_def_var(ncid, "east", NC_REAL, 1, &east_dimid, &east_varid));

    NC_CHECK(nc_def_var(ncid, "src_va", NC_REAL, 1, &depsrc_dimid, &src_va_varid));
    NC_CHECK(nc_def_var(ncid, "src_vb", NC_REAL, 1, &depsrc_dimid, &src_vb_varid));
    NC_CHECK(nc_def_var(ncid, "src_rho", NC_REAL, 1, &depsrc_dimid, &src_rho_varid));
    NC_CHECK(nc_def_var(ncid, "rcv_va", NC_REAL, 1, &deprcv_dimid, &rcv_va_varid));
    NC_CHECK(nc_def_var(ncid, "rcv_vb", NC_REAL, 1, &deprcv_dimid, &rcv_vb_varid));
    NC_CHECK(nc_def_var(ncid, "rcv_rho", NC_REAL, 1, &deprcv_dimid, &rcv_rho_varid));

    GRT_LOOP_ChnlGrid(im, c){
        int modr = GRT_SRC_M_ORDERS[im];
        char *s_title = NULL;
        if(modr == 0 && GRT_ZRT_CODES[c] == 'T') continue;

        GRT_SAFE_ASPRINTF(&s_title, "%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
        NC_CHECK(nc_def_var(ncid, s_title, NC_REAL, 4, dimids, &u_varids[im][c]));
        if(lib->calc_upar){
            GRT_SAFE_ASPRINTF(&s_title, "z%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
            NC_CHECK(nc_def_var(ncid, s_title, NC_REAL, 4, dimids, &uiz_varids[im][c]));
            GRT_SAFE_ASPRINTF(&s_title, "r%s%c", GRT_SRC_M_NAME_ABBR[im], GRT_ZRT_CODES[c]);
            NC_CHECK(nc_def_var(ncid, s_title, NC_REAL, 4, dimids, &uir_varids[im][c]));
        }
        GRT_SAFE_FREE_PTR(s_title);
    }

    NC_CHECK(nc_enddef(ncid));

    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, depsrc_varid, lib->depsrcs));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, deprcv_varid, lib->deprcvs));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, north_varid, lib->norths));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, east_varid, lib->easts));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, src_va_varid, lib->src_va));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, src_vb_varid, lib->src_vb));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, src_rho_varid, lib->src_rho));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, rcv_va_varid, lib->rcv_va));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, rcv_vb_varid, lib->rcv_vb));
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, rcv_rho_varid, lib->rcv_rho));

    size_t nr = lib->nnorth * lib->neast;
    real_t *tmpdata = (real_t *)calloc(nr, sizeof(real_t));
    for(size_t is = 0; is < lib->ndepsrc; ++is){
        for(size_t ir = 0; ir < lib->ndeprcv; ++ir){
            size_t start[4] = {is, ir, 0, 0};
            size_t count[4] = {1, 1, lib->nnorth, lib->neast};

            GRT_LOOP_ChnlGrid(im, c){
                int modr = GRT_SRC_M_ORDERS[im];
                if(modr == 0 && GRT_ZRT_CODES[c] == 'T') continue;

                for(size_t ipt = 0; ipt < nr; ++ipt){
                    tmpdata[ipt] = lib->u[is][ir][ipt][im][c];
                }
                NC_CHECK(NC_FUNC_REAL(nc_put_vara)(ncid, u_varids[im][c], start, count, tmpdata));

                if(lib->calc_upar){
                    for(size_t ipt = 0; ipt < nr; ++ipt){
                        tmpdata[ipt] = lib->uiz[is][ir][ipt][im][c];
                    }
                    NC_CHECK(NC_FUNC_REAL(nc_put_vara)(ncid, uiz_varids[im][c], start, count, tmpdata));
                    for(size_t ipt = 0; ipt < nr; ++ipt){
                        tmpdata[ipt] = lib->uir[is][ir][ipt][im][c];
                    }
                    NC_CHECK(NC_FUNC_REAL(nc_put_vara)(ncid, uir_varids[im][c], start, count, tmpdata));
                }
            }
        }
    }
    GRT_SAFE_FREE_PTR(tmpdata);
    NC_CHECK(nc_close(ncid));
}
