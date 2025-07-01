/**
 * @file   colorstr.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12
 * 
 * 定义一些字符串宏，以打印有颜色的字符串
 *                   
 */


#pragma once 

#define DEFAULT_RESTORE  "\033[0m" 
#define REGULAR_RED      "\033[0;31m"
#define BOLD_RED         "\033[1;31m"
#define REGULAR_GREEN    "\033[0;32m"
#define BOLD_GREEN       "\033[1;32m"
#define REGULAR_YELLOW   "\033[0;33m"
#define BOLD_YELLOW      "\033[1;33m"
#define REGULAR_BLUE     "\033[0;34m"
#define BOLD_BLUE        "\033[1;34m"
#define REGULAR_MAGENTA  "\033[0;35m"
#define BOLD_MAGENTA     "\033[1;35m"
#define REGULAR_CYAN     "\033[0;36m"
#define BOLD_CYAN        "\033[1;36m"