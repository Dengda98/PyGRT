/**
 * @file   stgrnlib.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-08
 *
 * 静态格林函数库 STGRNLIB：内存管理与四维 nc 读写
 *
 */

#include <string.h>
#include <math.h>

#include "grt/static/stgrnlib.h"
#include "grt/common/mynetcdf.h"
#include "grt/common/checkerror.h"
#include "grt/common/search.h"

/** 检查数组严格升序 */
static void require_strictly_ascending(const real_t *a, size_t n, const char *name)
{
    for(size_t i = 1; i < n; ++i){
        if(!(a[i] > a[i - 1])){
            GRTRaiseError("%s must be strictly ascending.", name);
        }
    }
}

/**
 * 由 norths/easts 填充网格序 rs，以及升序 sort_rs / sort_rs_idx / isUniform / dr
 *
 * isUniform 判定与 syn 一致：nr>2 且网格序已是等距升序
 */
static void fill_rs_meta(STGRNLIB *lib)
{
    lib->nr = lib->nnorth * lib->neast;
    lib->rs = (real_t *)malloc(lib->nr * sizeof(real_t));
    lib->sort_rs = (real_t *)malloc(lib->nr * sizeof(real_t));
    lib->sort_rs_idx = (size_t *)malloc(lib->nr * sizeof(size_t));

    for(size_t inorth = 0; inorth < lib->nnorth; ++inorth){
        for(size_t ieast = 0; ieast < lib->neast; ++ieast){
            size_t ipt = ieast + inorth * lib->neast;
            lib->rs[ipt] = hypot(lib->norths[inorth], lib->easts[ieast]);
            lib->sort_rs_idx[ipt] = ipt;
        }
    }
    memcpy(lib->sort_rs, lib->rs, lib->nr * sizeof(real_t));

    // 还未排序前，先判断是否是一个等距升序数组，便于加快后续查找
    lib->isUniform = (lib->nr > 2);
    lib->dr = (lib->nr > 1) ? (lib->rs[1] - lib->rs[0]) : 0.0;
    for(size_t ir = 1; ir + 1 < lib->nr; ++ir){
        if(fabs(2.0 * lib->rs[ir] - (lib->rs[ir - 1] + lib->rs[ir + 1])) > 1e-3
           || lib->rs[ir - 1] >= lib->rs[ir]
           || lib->rs[ir] >= lib->rs[ir + 1]){
            lib->isUniform = false;
            break;
        }
    }

    if(!lib->isUniform){
        if(lib->nr > 1
           && grt_argsort(lib->rs, lib->nr, sizeof(*lib->rs), grt_compare_real_t, lib->sort_rs_idx) != 0){
            GRTRaiseError("Unable to sort epicentral distances.");
        }
        for(size_t i = 0; i < lib->nr; ++i){
            lib->sort_rs[i] = lib->rs[lib->sort_rs_idx[i]];
        }
    }
}

/** 按已填好的维度申请 u/uiz/uir */
static void allocate_u(STGRNLIB *lib)
{
    size_t nr = lib->nr;
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
    require_strictly_ascending(depsrcs, ndepsrc, "depsrcs");
    require_strictly_ascending(deprcvs, ndeprcv, "deprcvs");
    require_strictly_ascending(norths, nnorth, "norths");
    require_strictly_ascending(easts, neast, "easts");

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

    fill_rs_meta(lib);

    // 介质参数由调用方随后填入
    lib->src_va = (real_t *)calloc(ndepsrc, sizeof(real_t));
    lib->src_vb = (real_t *)calloc(ndepsrc, sizeof(real_t));
    lib->src_rho = (real_t *)calloc(ndepsrc, sizeof(real_t));
    lib->rcv_va = (real_t *)calloc(ndeprcv, sizeof(real_t));
    lib->rcv_vb = (real_t *)calloc(ndeprcv, sizeof(real_t));
    lib->rcv_rho = (real_t *)calloc(ndeprcv, sizeof(real_t));
    lib->nlayer = 0;
    lib->modarr = NULL;

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
    GRT_SAFE_FREE_PTR(lib->rs);
    GRT_SAFE_FREE_PTR(lib->sort_rs);
    GRT_SAFE_FREE_PTR(lib->sort_rs_idx);
    GRT_SAFE_FREE_PTR(lib->src_va);
    GRT_SAFE_FREE_PTR(lib->src_vb);
    GRT_SAFE_FREE_PTR(lib->src_rho);
    GRT_SAFE_FREE_PTR(lib->rcv_va);
    GRT_SAFE_FREE_PTR(lib->rcv_vb);
    GRT_SAFE_FREE_PTR(lib->rcv_rho);
    GRT_SAFE_FREE_PTR(lib->modarr);
    free(lib);
}

void grt_stgrnlib_set_modarr(
    STGRNLIB *lib, size_t nlayer, const real_t (*modarr)[GRT_MODARR_NCOL])
{
    if(lib == NULL){
        GRTRaiseError("lib is NULL.");
    }
    if(nlayer == 0 || modarr == NULL){
        GRTRaiseError("nlayer and modarr must be non-empty.");
    }
    GRT_SAFE_FREE_PTR(lib->modarr);
    lib->modarr = (real_t (*)[GRT_MODARR_NCOL])malloc(
        sizeof(real_t) * GRT_MODARR_NCOL * nlayer);
    if(lib->modarr == NULL){
        GRTRaiseError("Failed to allocate modarr.");
    }
    memcpy(lib->modarr, modarr, sizeof(real_t) * GRT_MODARR_NCOL * nlayer);
    lib->nlayer = nlayer;
}

real_t grt_stgrnlib_default_subfault_size(const STGRNLIB *lib)
{
    const real_t atol = 1e-8;

    if(lib == NULL){
        GRTRaiseError("lib is NULL.");
    }
    if(lib->sort_rs == NULL){
        GRTRaiseError("STGRNLIB distance metadata is missing.");
    }

    real_t dr = -1.0;
    for(size_t i = 0; i + 1 < lib->nr; ++i){
        real_t d = lib->sort_rs[i + 1] - lib->sort_rs[i];
        if(d > atol && (dr < 0.0 || d < dr)) dr = d;
    }

    real_t dz = -1.0;
    for(size_t i = 0; i + 1 < lib->ndepsrc; ++i){
        real_t d = fabs(lib->depsrcs[i + 1] - lib->depsrcs[i]);
        if(d > atol && (dz < 0.0 || d < dz)) dz = d;
    }

    if(dr < 0.0 && dz < 0.0){
        GRTRaiseError(
            "Cannot infer default subfault size: "
            "need at least two distinct epicentral-distance or source-depth samples.");
    }
    if(dr < 0.0) return dz;
    if(dz < 0.0) return dr;
    return GRT_MIN(dr, dz);
}

/** 读入一层 (is, ir) 的通道；ncid 已打开，变量为 4D */
static void read_nc_channels_slice(STGRNLIB *lib, size_t is, size_t ir, int ncid)
{
    size_t nr = lib->nr;
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

    // 建库模型矩阵（必需）
    {
        int layer_dimid, param_dimid;
        size_t nlayer = 0, nparam = 0;
        if(nc_inq_dimid(ncid, "layer", &layer_dimid) != NC_NOERR ||
           nc_inq_dimid(ncid, "model_param", &param_dimid) != NC_NOERR){
            NC_CHECK(nc_close(ncid));
            grt_stgrnlib_free(lib);
            GRTRaiseError(
                "Invalid STGRNLIB nc \"%s\": missing model dimensions "
                "(layer, model_param); rebuild with current greenfn.",
                path);
        }
        NC_CHECK(nc_inq_dimlen(ncid, layer_dimid, &nlayer));
        NC_CHECK(nc_inq_dimlen(ncid, param_dimid, &nparam));
        if(nlayer == 0 || nparam != GRT_MODARR_NCOL){
            NC_CHECK(nc_close(ncid));
            grt_stgrnlib_free(lib);
            GRTRaiseError(
                "Invalid STGRNLIB nc \"%s\": layer=%zu, model_param=%zu "
                "(expect layer>0, model_param=%d).",
                path, nlayer, nparam, GRT_MODARR_NCOL);
        }
        NC_CHECK(nc_inq_varid(ncid, "model", &varid));
        real_t (*modarr)[GRT_MODARR_NCOL] = (real_t (*)[GRT_MODARR_NCOL])malloc(
            sizeof(real_t) * GRT_MODARR_NCOL * nlayer);
        NC_CHECK(NC_FUNC_REAL(nc_get_var)(ncid, varid, (real_t *)modarr));
        grt_stgrnlib_set_modarr(lib, nlayer, (const real_t (*)[GRT_MODARR_NCOL])modarr);
        GRT_SAFE_FREE_PTR(modarr);
    }

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
    if(lib->nlayer == 0 || lib->modarr == NULL){
        GRTRaiseError("STGRNLIB has no model matrix; call grt_stgrnlib_set_modarr first.");
    }

    int ncid;
    int depsrc_dimid, deprcv_dimid, north_dimid, east_dimid;
    int layer_dimid, param_dimid;
    int dimids[4];
    int depsrc_varid, deprcv_varid, north_varid, east_varid;
    int model_varid;
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
    NC_CHECK(nc_def_dim(ncid, "layer", lib->nlayer, &layer_dimid));
    NC_CHECK(nc_def_dim(ncid, "model_param", GRT_MODARR_NCOL, &param_dimid));
    dimids[0] = depsrc_dimid;
    dimids[1] = deprcv_dimid;
    dimids[2] = north_dimid;
    dimids[3] = east_dimid;

    {
        int model_dimids[2] = {layer_dimid, param_dimid};
        NC_CHECK(nc_def_var(ncid, "model", NC_REAL, 2, model_dimids, &model_varid));
    }

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
    NC_CHECK(NC_FUNC_REAL(nc_put_var)(ncid, model_varid, (const real_t *)lib->modarr));

    size_t nr = lib->nr;
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
