/**
 * @file   util.c
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2025-08
 * 
 * 其它辅助函数
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

#include "grt/common/util.h"
#include "grt/common/const.h"
#include "grt/common/checkerror.h"

char ** grt_string_split(const char *string, const char *delim, size_t *size)
{
    char *str_copy = strdup(string);  // 创建字符串副本，以免修改原始字符串
    char *token = strtok(str_copy, delim);

    char **s_split = NULL;
    *size = 0;

    while(token != NULL){
        s_split = (char**)realloc(s_split, sizeof(char*)*(*size+1));
        s_split[*size] = NULL;
        s_split[*size] = (char*)realloc(s_split[*size], sizeof(char)*(strlen(token)+1));
        strcpy(s_split[*size], token);

        token = strtok(NULL, delim);
        (*size)++;
    }
    free(str_copy);

    return s_split;
}

char ** grt_string_from_file(FILE *fp, size_t *size){
    char **s_split = NULL;
    *size = 0;
    s_split = (char**)realloc(s_split, sizeof(char*)*(*size+1));
    s_split[*size] = NULL;

    size_t len=0;
    while(grt_getline(&s_split[*size], &len, fp) != -1){
        size_t line_len = strlen(s_split[*size]);
        if(line_len > 0 && s_split[*size][line_len - 1] == '\n'){
            s_split[*size][--line_len] = '\0';
        }
        if(line_len > 0 && s_split[*size][line_len - 1] == '\r'){
            s_split[*size][line_len - 1] = '\0';
        }
        (*size)++;
        s_split = (char**)realloc(s_split, sizeof(char*)*(*size+1));
        s_split[*size] = NULL;
    }
    return s_split;
}

real_t *grt_parse_real_array(const char *optarg, size_t *size, char ***s_values, char optname)
{
    real_t a1, a2, delta;
    char **s_vals = NULL;
    size_t n = 0;

    // 逗号分隔的数值列表
    if(grt_string_composed_of(optarg, GRT_NUM_STR "eE+-" ".,")){
        s_vals = grt_string_split(optarg, ",", &n);
    }
    // 等间距范围
    else if(3 == sscanf(optarg, "%lf/%lf/%lf", &a1, &a2, &delta)){
        if(delta <= 0.0){
            GRTRaiseError("-%c: nonpositive spacing (%f).", optname, delta);
        }
        if(a1 > a2){
            GRTRaiseError("-%c: start (%f) > end (%f).", optname, a1, a2);
        }
        n = (size_t)floor((a2 - a1) / delta) + 1;
        s_vals = (char **)calloc(n, sizeof(char *));
        for(size_t i = 0; i < n; ++i){
            GRT_SAFE_ASPRINTF(&s_vals[i], "%.15g", a1 + delta * i);
        }
    }
    // 文件中的每行读取一个数值
    else {
        FILE *fp = GRTCheckOpenFile(optarg, "r");
        s_vals = grt_string_from_file(fp, &n);
        fclose(fp);
    }

    if(n == 0){
        GRTRaiseError("-%c: empty value list.", optname);
    }

    real_t *values = (real_t *)calloc(n, sizeof(real_t));
    for(size_t i = 0; i < n; ++i){
        values[i] = atof(s_vals[i]);
        if(values[i] < 0.0){
            GRTRaiseError("-%c: negative value (%f) is not supported.", optname, values[i]);
        }
        if(i > 0 && !(values[i] > values[i - 1])){
            GRTRaiseError("-%c: values must be strictly ascending.", optname);
        }
    }

    if(s_values != NULL){
        *s_values = s_vals;
    } else {
        GRT_SAFE_FREE_PTR_ARRAY(s_vals, n);
    }
    *size = n;
    return values;
}

bool grt_string_composed_of(const char *str, const char *alws){
    bool allowed[256] = {false};  // 初始全为false（不允许）

    // 标记允许的字符
    for (int i = 0; alws[i] != '\0'; i++) {
        unsigned char c = alws[i];  // 转为无符号避免负数索引
        allowed[c] = true;
    }

    // 检查目标字符串中的每个字符
    for (int i = 0; str[i] != '\0'; i++) {
        unsigned char c = str[i];
        if (!allowed[c]) {  // 若字符不在允许集合中
            return false;
        }
    }

    // 所有字符均在允许集合中
    return true;
}

int grt_string_ncols(const char *string, const char* delim){
    int count = 0;
    
    const char *str = string;
    while (*str) {
        // 跳过所有分隔符
        while (*str && strchr(delim, *str)) str++;
        // 如果还有非分隔符字符，增加计数
        if (*str) count++;
        // 跳过所有非分隔符字符
        while (*str && !strchr(delim, *str)) str++;
    }
    
    return count;
}


const char* grt_get_basename(const char* path) {
    // 找到最后一个 '/'
    char* last_slash = strrchr(path, '/'); 
    
#ifdef _WIN32
    char* last_backslash = strrchr(path, '\\');
    if (last_backslash && (!last_slash || last_backslash > last_slash)) {
        last_slash = last_backslash;
    }
#endif
    if (last_slash) {
        // 返回最后一个 '/' 之后的部分
        return last_slash + 1; 
    }
    // 如果没有 '/'，整个路径就是最后一项
    return path; 
}


/**
 * 判断指定字符串区间是否表示一个有限实数
 *
 * 字符串区间左闭右开，调用方无需为区间补充字符串结束符
 *
 * @param[in]    begin    字符串区间的起始位置，包含该位置
 * @param[in]    end      字符串区间的结束位置，不包含该位置
 *
 * @return   如果区间是一个有限实数则返回 true，否则返回 false
 */
static bool grt_is_numeric_field(const char *begin, const char *end)
{
    // 忽略数值字段两端的空白字符
    while(begin < end && isspace((unsigned char)*begin)) ++begin;
    while(end > begin && isspace((unsigned char)end[-1])) --end;
    if(begin == end) return false;

    // strtod 需要以字符串结束符结尾，因此复制指定的字符串区间
    size_t length = (size_t)(end - begin);
    char *value = (char *)calloc(length + 1, sizeof(char));
    memcpy(value, begin, length);

    // strtod 支持小数、科学计数法以及正负号
    errno = 0;
    char *value_end = NULL;
    real_t number = strtod(value, &value_end);

    // 要求完整解析字段，并排除溢出、下溢、无穷大和非数值
    bool valid = (value_end != value) && (*value_end == '\0')
        && (errno != ERANGE) && isfinite(number);

    GRT_SAFE_FREE_PTR(value);
    return valid;
}


/**
 * 判断格林函数子目录名是否符合当前模型的命名格式
 *
 * 子目录名格式为 modelname_depsrc_deprcv_dist
 * 数值字段从右侧分隔，避免模型名中的下划线参与解析
 *
 * @param[in]    name       待检查的子目录名
 * @param[in]    modelname  当前模型文件名
 *
 * @return   如果名称符合当前模型的格林函数目录格式则返回 true，否则返回 false
 */
static bool grt_is_greenfn_subdir_name(const char *name, const char *modelname)
{
    size_t model_length = strlen(modelname);
    size_t name_length = strlen(name);

    // 先精确匹配模型名，并要求模型名后紧跟一个下划线
    if(name_length <= model_length || strncmp(name, modelname, model_length) != 0 || name[model_length] != '_'){
        return false;
    }

    const char *suffix = name + model_length + 1;
    size_t suffix_length = name_length - model_length - 1;

    // 从右侧寻找震中距和接收器深度之间的分隔符
    const char *separator2_ptr = strrchr(suffix, '_');
    if(separator2_ptr == NULL) return false;

    size_t separator2 = (size_t)(separator2_ptr - suffix);
    size_t separator1 = separator2;
    // 继续向左寻找震源深度和接收器深度之间的分隔符
    while(separator1 > 0 && suffix[separator1 - 1] != '_') --separator1;
    if(separator1 == 0) return false;
    --separator1;

    return grt_is_numeric_field(suffix, suffix + separator1)
        && grt_is_numeric_field(suffix + separator1 + 1, suffix + separator2)
        && grt_is_numeric_field(suffix + separator2 + 1, suffix + suffix_length);
}


/**
 * 检查动态和模态格林函数输出根目录中的文件和目录
 *
 * 根目录允许包含 command、当前模型副本，以及当前模型对应的格林函数子目录
 * 调用方应在执行本函数前创建 output_dir
 *
 * @param[in]    output_dir  格林函数输出根目录
 * @param[in]    modelname   当前模型文件名
 */
void grt_check_greenfn_output_dir(const char *output_dir, const char *modelname)
{
    // 目录应已由调用方创建，这里打开目录并检查其中的已有条目
    DIR *dir = opendir(output_dir);
    if(dir == NULL){
        GRTRaiseError("Cannot open Green's function output directory \"%s\". Error code: %d\n", output_dir, errno);
    }

    errno = 0;
    struct dirent *entry;
    while((entry = readdir(dir)) != NULL){
        // 跳过目录自身和父目录项
        if(strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0){
            continue;
        }

        // 使用完整路径检查条目类型
        char *entry_path = NULL;
        GRT_SAFE_ASPRINTF(&entry_path, "%s/%s", output_dir, entry->d_name);

        struct stat entry_stat;
        if(stat(entry_path, &entry_stat) != 0){
            int error_code = errno;
            closedir(dir);
            GRTRaiseError("Cannot inspect Green's function output entry \"%s\". Error code: %d\n", entry_path, error_code);
        }

        // 根目录只允许 command、模型副本和当前模型的格林函数子目录
        bool allowed = false;
        if(strcmp(entry->d_name, "command") == 0 || strcmp(entry->d_name, modelname) == 0){
            allowed = S_ISREG(entry_stat.st_mode);
        } else if(S_ISDIR(entry_stat.st_mode)){
            allowed = grt_is_greenfn_subdir_name(entry->d_name, modelname);
        }

        if(!allowed){
            closedir(dir);
            GRTRaiseError("The current model is \"%s\", but the output directory contains \"%s\". This is not allowed.\n", modelname, entry_path);
        }

        GRT_SAFE_FREE_PTR(entry_path);
    }

    // readdir 返回 NULL 时，通过 errno 区分正常结束和读取失败
    if(errno != 0){
        int error_code = errno;
        closedir(dir);
        GRTRaiseError("Cannot read Green's function output directory \"%s\". Error code: %d\n", output_dir, error_code);
    }

    closedir(dir);
}


void grt_trim_whitespace(char* str) {
    char* end;
    
    // 去除首部空白
    while (isspace((unsigned char)*str)) str++;
    
    // 去除尾部空白
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    
    // 写入终止符
    *(end + 1) = '\0';
}


bool grt_is_empty_line(const char* line) {
    // 跳过前导空白
    while (isspace((unsigned char)*line)) line++;

    // 检查是否为空行
    return *line == '\0';
}


bool grt_is_comment_line(const char* line) {
    // 跳过前导空白
    while (isspace((unsigned char)*line)) line++;

    // 检查是否为注释行
    return *line == GRT_COMMENT_HEAD;
}


bool grt_is_comment_or_empty_line(const char* line) {
    return grt_is_empty_line(line) || grt_is_comment_line(line);
}


void grt_copy_file(const char *src, const char *dst)
{
    // 源与目标为同一路径时跳过，避免截断自身
    {
#if _TEST_WHETHER_WIN32_
        char src_full[_MAX_PATH];
        char dst_full[_MAX_PATH];
        if(_fullpath(src_full, src, _MAX_PATH) != NULL
           && _fullpath(dst_full, dst, _MAX_PATH) != NULL
           && strcmp(src_full, dst_full) == 0){
            return;
        }
#else
        char src_real[PATH_MAX];
        char dst_real[PATH_MAX];
        if(realpath(src, src_real) != NULL && realpath(dst, dst_real) != NULL
           && strcmp(src_real, dst_real) == 0){
            return;
        }
#endif
    }

    FILE *fin = GRTCheckOpenFile(src, "rb");
    FILE *fout = GRTCheckOpenFile(dst, "wb");

    char buf[8192];
    size_t n;
    while((n = fread(buf, 1, sizeof(buf), fin)) > 0){
        if(fwrite(buf, 1, n, fout) != n){
            fclose(fin);
            fclose(fout);
            GRTRaiseError("Failed to write copy of \"%s\" to \"%s\".\n", src, dst);
        }
    }
    if(ferror(fin)){
        fclose(fin);
        fclose(fout);
        GRTRaiseError("Failed to read \"%s\" while copying to \"%s\".\n", src, dst);
    }

    fclose(fin);
    fclose(fout);
}


ssize_t grt_getline(char **lineptr, size_t *n, FILE *stream){
    if (!lineptr || !n || !stream) {
        return -1;
    }
    
    char *buf = *lineptr;
    size_t size = *n;
    size_t len = 0;
    int c;
    
    // 如果缓冲区为空，分配初始缓冲区
    if (buf == NULL || size == 0) {
        size = 128;
        buf = malloc(size);
        if (buf == NULL) {
            return -1;
        }
    }
    
    // 逐字符读取直到换行符或EOF
    while ((c = fgetc(stream)) != EOF) {
        // 检查是否需要扩展缓冲区
        if (len + 1 >= size) {
            size_t new_size = size * 2;
            char *new_buf = realloc(buf, new_size);
            if (new_buf == NULL) {
                free(buf);
                return -1;
            }
            buf = new_buf;
            size = new_size;
        }
        
        buf[len++] = c;
        
        // 遇到换行符停止读取
        if (c == '\n') {
            break;
        }
    }
    
    // 如果没有读取到任何字符且遇到EOF
    if (len == 0 && c == EOF) {
        return -1;
    }
    
    // 添加字符串终止符
    buf[len] = '\0';
    
    *lineptr = buf;
    *n = size;
    
    return len;
}
