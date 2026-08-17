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


bool grt_is_comment_or_empty(const char* line) {
    // 跳过前导空白
    while (isspace((unsigned char)*line)) line++;
    
    // 检查是否为空行或注释行
    return (*line == '\0' || *line == GRT_COMMENT_HEAD);
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
