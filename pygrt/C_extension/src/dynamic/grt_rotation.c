/**
 * @file   grt_rotation.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-04-21
 * 
 *    根据预先合成的位移空间导数，组合成旋转张量
 * 
 */


#include "grt/common/sacio2.h"
#include "grt/common/const.h"

#include "grt.h"

/** 该子模块的参数控制结构体 */
typedef struct {
    char *name;
    char *s_synpath;
} GRT_MODULE_CTRL;


/** 释放结构体的内存 */
static void free_Ctrl(GRT_MODULE_CTRL *Ctrl){
    GRT_SAFE_FREE_PTR(Ctrl->name);
    GRT_SAFE_FREE_PTR(Ctrl->s_synpath);
    GRT_SAFE_FREE_PTR(Ctrl);
}


/** 打印使用说明 */
static void print_help(){
printf("\n"
"[grt rotation] %s\n\n", GRT_VERSION);printf(
"    Conbine spatial derivatives of displacements into rotation tensor.\n"
"    For example, \"ZR\" in filename means 0.5*(u_{z,r} - u_{r,z}).\n"
"\n\n"
"Usage:\n"
"----------------------------------------------------------------\n"
"    grt rotation <syn_dir>\n"
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
int rotation_main(int argc, char **argv){
    GRT_MODULE_CTRL *Ctrl = calloc(1, sizeof(*Ctrl));
    Ctrl->name = strdup(argv[0]);
    const char *command = Ctrl->name;

    getopt_from_command(Ctrl, argc, argv);
    
    // 合成地震图目录路径
    Ctrl->s_synpath = strdup(argv[1]);

    // 检查是否存在该目录
    GRTCheckDirExist(command, Ctrl->s_synpath);

    // ----------------------------------------------------------------------------------
    // 开始读取计算，输出3个量
    float *arrin = NULL;
    char c1, c2;
    char *s_filepath = NULL;

    // 输出分量格式，即是否需要旋转到ZNE
    bool rot2ZNE = false;
    // 三分量
    const char *chs = NULL;

    // 判断标志性文件是否存在，来判断输出使用ZNE还是ZRT
    GRT_SAFE_ASPRINTF(&s_filepath, "%s/nN.sac", Ctrl->s_synpath);
    rot2ZNE = (access(s_filepath, F_OK) == 0);

    // 指示特定的通道名
    chs = (rot2ZNE)? GRT_ZNE_CODES : GRT_ZRT_CODES;


    // 读取一个头段变量，获得基本参数，分配数组内存
    SACHEAD hd;
    GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->s_synpath, tolower(chs[0]), chs[0]);
    grt_read_SAC_HEAD(command, s_filepath, &hd);
    int npts=hd.npts;
    float dist=hd.dist;
    float *arrout = (float*)calloc(npts, sizeof(float));

    // ----------------------------------------------------------------------------------
    // 循环3个分量
    for(int i1=0; i1<2; ++i1){
        c1 = chs[i1];
        for(int i2=i1+1; i2<3; ++i2){
            c2 = chs[i2];

            // 读取数据 u_{i,j}
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->s_synpath, tolower(c2), c1);
            arrin = grt_read_SAC(command, s_filepath, &hd, arrin);

            // 累加
            for(int i=0; i<npts; ++i)  arrout[i] += arrin[i];

            // 读取数据 u_{j,i}
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/%c%c.sac", Ctrl->s_synpath, tolower(c1), c2);
            arrin = grt_read_SAC(command, s_filepath, &hd, arrin);

            // 累加
            for(int i=0; i<npts; ++i)  arrout[i] = (arrout[i] - arrin[i]) * 0.5f;

            // 特殊情况需加上协变导数，1e-5是因为km->cm
            if(c1=='R' && c2=='T'){
                // 读取数据 u_T
                GRT_SAFE_ASPRINTF(&s_filepath, "%s/T.sac", Ctrl->s_synpath);
                arrin = grt_read_SAC(command, s_filepath, &hd, arrin);
                for(int i=0; i<npts; ++i)  arrout[i] -= 0.5f * arrin[i] / dist * 1e-5;
            }

            // 保存到SAC
            sprintf(hd.kcmpnm, "%c%c", c1, c2);
            GRT_SAFE_ASPRINTF(&s_filepath, "%s/rotation_%c%c.sac", Ctrl->s_synpath, c1, c2);
            write_sac(s_filepath, hd, arrout);

            // 置零
            for(int i=0; i<npts; ++i)  arrout[i] = 0.0f;
        }
    }

    GRT_SAFE_FREE_PTR(arrin);
    GRT_SAFE_FREE_PTR(arrout);
    GRT_SAFE_FREE_PTR(s_filepath);

    free_Ctrl(Ctrl);
    return EXIT_SUCCESS;
}
