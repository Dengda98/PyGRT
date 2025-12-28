/**
 * @file   grt_travt.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-08
 * 
 *    计算一维均匀半无限层状介质的初至走时的主函数
 * 
 */

#include "grt/common/travt.h"
#include "grt/common/const.h"
#include "grt/common/model.h"
#include "grt/common/util.h"

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    /** 输入模型 */
    struct {
        bool active;
        char *s_modelpath;
        GRT_MODEL1D *mod1d;         ///< 模型结构体指针
    } M;
    /** 震源和接收器深度 */
    struct {
        bool active;
        real_t depsrc;
        real_t deprcv;
        char *s_depsrc;
        char *s_deprcv;
    } D;
    /** 震中距 */
    struct {
        bool active;
        char **s_rs;
        real_t *rs;
        size_t nr;
    } R;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->M.s_modelpath);
    grt_free_mod1d(Ctrl->M.mod1d);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_depsrc);
    GRT_SAFE_FREE_PTR(Ctrl->D.s_deprcv);

    GRT_SAFE_FREE_PTR(Ctrl->R.rs);
    GRT_SAFE_FREE_PTR_ARRAY(Ctrl->R.s_rs, Ctrl->R.nr);

    GRT_SAFE_FREE_PTR(Ctrl);
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
                    GRTBadOptionError(D, "");
                };
                if(1 != sscanf(Ctrl->D.s_depsrc, "%lf", &Ctrl->D.depsrc)){
                    GRTBadOptionError(D, "");
                }
                if(1 != sscanf(Ctrl->D.s_deprcv, "%lf", &Ctrl->D.deprcv)){
                    GRTBadOptionError(D, "");
                }
                if(Ctrl->D.depsrc < 0.0 || Ctrl->D.deprcv < 0.0){
                    GRTBadOptionError(D, "Negative value in -D is not supported.");
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
                    FILE *fp = GRTCheckOpenFile(optarg, "r");
                    Ctrl->R.s_rs = grt_string_from_file(fp, &Ctrl->R.nr);
                    fclose(fp);
                }
                // 转为浮点数
                Ctrl->R.rs = (real_t*)realloc(Ctrl->R.rs, sizeof(real_t)*(Ctrl->R.nr));
                for(size_t i=0; i<Ctrl->R.nr; ++i){
                    Ctrl->R.rs[i] = atof(Ctrl->R.s_rs[i]);
                    if(Ctrl->R.rs[i] < 0.0){
                        GRTBadOptionError(R, "Can't set negative epicentral distance(%f).", Ctrl->R.rs[i]);
                    }
                }
                break;

            GRT_Common_Options_in_Switch((char)(optopt));
        }
    }

    // 检查必须设置的参数是否有设置
    GRTCheckOptionSet(argc > 1);
    GRTCheckOptionActive(Ctrl, M);
    GRTCheckOptionActive(Ctrl, D);
    GRTCheckOptionActive(Ctrl, R);

}


/** 子模块主函数 */
int travt_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));

    getopt_from_command(Ctrl, argc, argv);

    // 读入模型文件
    if((Ctrl->M.mod1d = grt_read_mod1d_from_file(Ctrl->M.s_modelpath, Ctrl->D.depsrc, Ctrl->D.deprcv, true)) == NULL){
        exit(EXIT_FAILURE);
    }
    GRT_MODEL1D *mod1d = Ctrl->M.mod1d;

    printf("------------------------------------------------\n");
    printf(" Distance(km)     Tp(secs)         Ts(secs)     \n");
    real_t travtP=-1, travtS=-1;
    for(size_t i=0; i<Ctrl->R.nr; ++i){
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