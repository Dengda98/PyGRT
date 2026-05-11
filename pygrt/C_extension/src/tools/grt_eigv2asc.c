/**
 * @file   grt_eigv2asc.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2026-02
 * 
 *     将 eigenv 模块输出的频散结果转为文本格式以供阅读
 * 
 */

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入文件路径 */
    char *s_filepath;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->s_filepath);
    GRT_SAFE_FREE_PTR(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt eigv2asc] %s\n\n", GRT_VERSION);printf(
"    Convert a .nc file from `eigenv` module into an ASCII file, \n"
"    write to standard output (ignore some vars).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt eigv2asc <ncfile>\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    (void)Ctrl;
    int opt;
    while ((opt = getopt(argc, argv, ":h")) != -1) {
        switch (opt) {
            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(argc > 1);
}



int eigv2asc_main(int argc, char **argv)
{
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    Ctrl->s_filepath = strdup(argv[1]);
    
    // 检查文件名是否存在
    GRTCheckFileExist(Ctrl->s_filepath);

    int ncid;
    int f_dimid, n_dimid;
    int f_varid;
    int c_varid, ciref_varid, cnum_varid;

    // 打开 NC 文件
    GRTCheckFileExist(Ctrl->s_filepath);
    NC_CHECK(nc_open(Ctrl->s_filepath, NC_NOWRITE, &ncid));

    size_t nfreq, nmode;

    // 读取 NC 文件的坐标变量
    NC_CHECK(nc_inq_dimid(ncid, "freq", &f_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, f_dimid, &nfreq));
    NC_CHECK(nc_inq_dimid(ncid, "mode", &n_dimid));
    NC_CHECK(nc_inq_dimlen(ncid, n_dimid, &nmode));
    real_t *freqs = (real_t *)calloc(nfreq, sizeof(real_t));
    NC_CHECK(nc_inq_varid(ncid, "freq", &f_varid));
    NC_CHECK(NC_FUNC_REAL(nc_get_var) (ncid, f_varid, freqs));

    // 每个频率下有多少阶
    int *cnum = (int *)calloc(nfreq, sizeof(int));
    NC_CHECK(nc_inq_varid(ncid, "cnum", &cnum_varid));
    NC_CHECK(NC_FUNC_INT(nc_get_var) (ncid, cnum_varid, cnum));

    NC_CHECK(nc_inq_varid(ncid, "c", &c_varid));
    NC_CHECK(nc_inq_varid(ncid, "ciref", &ciref_varid));    

    // 读入所有频散
    const int ndims = 2;
    size_t startp[ndims];
    size_t countp[ndims];
    for(size_t iw = 0; iw < nfreq; ++iw){
        int num = cnum[iw];
        real_t *c_roots = (real_t *)calloc(cnum[iw], sizeof(real_t));
        int *c_roots_iref = (int *)calloc(cnum[iw], sizeof(int));

        startp[0] = iw;
        startp[1] = 0;
        countp[0] = 1;
        countp[1] = cnum[iw];

        NC_CHECK(NC_FUNC_REAL(nc_get_vara) (ncid, c_varid, startp, countp, c_roots));
        NC_CHECK(NC_FUNC_INT(nc_get_vara) (ncid, ciref_varid, startp, countp, c_roots_iref));

        // 打印
        for(int ik = 0; ik < num; ++ik){
            printf("%15.5e %20.8e %10d %10d\n", freqs[iw], c_roots[ik], ik, c_roots_iref[ik]);
        }

        GRT_SAFE_FREE_PTR(c_roots);
        GRT_SAFE_FREE_PTR(c_roots_iref);
    }

    // 关闭文件
    NC_CHECK(nc_close(ncid));

    GRT_SAFE_FREE_PTR(freqs);
    GRT_SAFE_FREE_PTR(cnum);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}