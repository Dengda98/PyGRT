/**
 * @file   logo.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12
 * 
 *    logo字符串
 * 
 */

#pragma once


#include <stdio.h>

#include "const.h"
#include "version.h"
#include "colorstr.h"


inline GCC_ALWAYS_INLINE void print_logo(){
printf(BOLD_GREEN "\n"
"╔═══════════════════════════════════════════════════════════════╗\n"
"║                                                               ║\n"
"║              ██████╗     ██████╗     ████████╗                ║\n"
"║             ██╔════╝     ██╔══██╗    ╚══██╔══╝                ║\n"
"║             ██║  ███╗    ██████╔╝       ██║                   ║\n"
"║             ██║   ██║    ██╔══██╗       ██║                   ║\n"
"║             ╚██████╔╝    ██║  ██║       ██║                   ║\n"
"║              ╚═════╝     ╚═╝  ╚═╝       ╚═╝                   ║\n"
"║                                                               ║\n"
"║                                                               ║\n"
"║               Author: Zhu Dengda                              ║\n"
"║                Email: zhudengda@mail.iggcas.ac.cn             ║\n"
"║        Code Homepage: https://github.com/Dengda98/PyGRT       ║\n"
"║              License: GPL-3.0 license                         ║\n"
"║              Version: %-20s                    ║\n"
"║                                                               ║\n"
"║                                                               ║\n"
"║    A Command-Line Tool for Computing Synthetic Seismograms    ║\n"
"║            in Horizontally Layered Halfspace Model,           ║\n"
"║     using Generalized Reflection-Transmission Method(GRTM)    ║\n"
"║              and Discrete Wavenumber Method(DWM).             ║\n"
"║                                                               ║\n"
"╚═══════════════════════════════════════════════════════════════╝\n" 
DEFAULT_RESTORE, GRT_VERSION);

}