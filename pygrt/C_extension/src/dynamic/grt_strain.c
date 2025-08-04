/**
 * @file   grt_strain.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-03-28
 * 
 *    根据预先合成的位移空间导数，组合成应变张量
 * 
 */


#include "common/sacio2.h"
#include "common/const.h"

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    char *name;
    char *s_dirpath;
    char *s_prefix;
    char *s_synpath;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    free(Ctrl->name);
    free(Ctrl->s_dirpath);
    free(Ctrl->s_prefix);
    free(Ctrl->s_synpath);
    free(Ctrl);
}

/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt strain] %s\n\n", GRT_VERSION);printf(
"    Conbine spatial derivatives of displacements into strain tensor.\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt strain <syn_dir>/<name>\n"
"\n\n\n"
);
}


/** 从命令行中读取选项，处理后记录到全局变量中 */
static void getopt_from_command(GRT_MODULE_CTRL *Ctrl, int argc, char **argv){
    char* command = Ctrl->name;
    int opt;
    while ((opt = getopt(argc, argv, ":h")) != -1) {
        switch (opt) {
            GRT_Common_Options_in_Switch(command, (char)(optopt));
        }
    }

    // 检查必选项有没有设置
    GRTCheckOptionSet(command, argc > 1);
}



/** 子模块主函数 */
int strain_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);
    Ctrl->s_dirpath = strdup(argv[1]);
    
    const char *command = Ctrl->name;

    getopt_from_command(Ctrl, argc, argv);

    
    // 合成地震图目录路径
    Ctrl->s_synpath = (char*)malloc(sizeof(char)*(strlen(Ctrl->s_dirpath)+1));
    // 保存文件前缀 
    Ctrl->s_prefix = (char*)malloc(sizeof(char)*(strlen(Ctrl->s_dirpath)+1));
    if(2 != sscanf(Ctrl->s_dirpath, "%[^/]/%s", Ctrl->s_synpath, Ctrl->s_prefix)){
        GRTRaiseError("[%s] " BOLD_RED "Error format in \"%s\".\n" DEFAULT_RESTORE, command, Ctrl->s_dirpath);
    }

    // 检查是否存在该目录
    GRTCheckDirExist(command, Ctrl->s_synpath);

    // ----------------------------------------------------------------------------------
    // 开始读取计算，输出6个量
    float *arrin = NULL;
    char c1, c2;
    char *s_filepath = (char*)malloc(sizeof(char) * (strlen(Ctrl->s_synpath)+strlen(Ctrl->s_prefix)+100));

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = false;
    // 三分量
    const char *chs = NULL;

    // 判断标志性文件是否存在，来判断输出使用ZNE还是ZRT
    sprintf(s_filepath, "%s/n%sN.sac", Ctrl->s_synpath, Ctrl->s_prefix);
    rot2ZNE = (access(s_filepath, F_OK) == 0);

    // 指示特定的通道名
    chs = (rot2ZNE)? ZNEchs : ZRTchs;


    // 读取一个头段变量，获得基本参数，分配数组内存
    SACHEAD hd;
    sprintf(s_filepath, "%s/%c%s%c.sac", Ctrl->s_synpath, tolower(chs[0]), Ctrl->s_prefix, chs[0]);
    read_SAC_HEAD(command, s_filepath, &hd);
    int npts=hd.npts;
    float dist=hd.dist;
    float *arrout = (float*)calloc(npts, sizeof(float));

    // ----------------------------------------------------------------------------------
    // 循环6个分量
    for(int i1=0; i1<3; ++i1){
        c1 = chs[i1];
        for(int i2=i1; i2<3; ++i2){
            c2 = chs[i2];

            // 读取数据 u_{i,j}
            sprintf(s_filepath, "%s/%c%s%c.sac", Ctrl->s_synpath, tolower(c2), Ctrl->s_prefix, c1);
            arrin = read_SAC(command, s_filepath, &hd, arrin);

            // 累加
            for(int i=0; i<npts; ++i)  arrout[i] += arrin[i];

            // 读取数据 u_{j,i}
            sprintf(s_filepath, "%s/%c%s%c.sac", Ctrl->s_synpath, tolower(c1), Ctrl->s_prefix, c2);
            arrin = read_SAC(command, s_filepath, &hd, arrin);

            // 累加
            for(int i=0; i<npts; ++i)  arrout[i] = (arrout[i] + arrin[i]) * 0.5f;

            // 特殊情况需加上协变导数，1e-5是因为km->cm
            if(c1=='R' && c2=='T'){
                // 读取数据 u_T
                sprintf(s_filepath, "%s/%sT.sac", Ctrl->s_synpath, Ctrl->s_prefix);
                arrin = read_SAC(command, s_filepath, &hd, arrin);
                for(int i=0; i<npts; ++i)  arrout[i] -= 0.5f * arrin[i] / dist * 1e-5;
            }
            else if(c1=='T' && c2=='T'){
                // 读取数据 u_R
                sprintf(s_filepath, "%s/%sR.sac", Ctrl->s_synpath, Ctrl->s_prefix);
                arrin = read_SAC(command, s_filepath, &hd, arrin);
                for(int i=0; i<npts; ++i)  arrout[i] += arrin[i] / dist * 1e-5;
            }

            // 保存到SAC
            sprintf(hd.kcmpnm, "%c%c", c1, c2);
            sprintf(s_filepath, "%s/%s.strain.%c%c.sac", Ctrl->s_synpath, Ctrl->s_prefix, c1, c2);
            write_sac(s_filepath, hd, arrout);

            // 置零
            for(int i=0; i<npts; ++i)  arrout[i] = 0.0f;
        }
    }

    if(arrin)   free(arrin);
    if(arrout)  free(arrout);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
