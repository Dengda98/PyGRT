/**
 * @file   matrix.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-07-24
 * 
 * 2x2小矩阵的加、减、乘、除、求逆等操作，由于均为小型数组操作，所有函数均为内联函数               
 */

#pragma once

#include "common/const.h"

/**
 * 计算2x2复矩阵的逆  
 * 
 * @param[in]      M        原矩阵
 * @param[out]     invM     逆矩阵
 * @param[out]     stats    状态代码，是否有除零错误，非0为异常值
 */ 
inline GCC_ALWAYS_INLINE void cmat2x2_inv(const MYCOMPLEX M[2][2], MYCOMPLEX invM[2][2], MYINT *stats) {
    MYCOMPLEX M00 = M[0][0];
    MYCOMPLEX M01 = M[0][1];
    MYCOMPLEX M10 = M[1][0];
    MYCOMPLEX M11 = M[1][1];
    MYCOMPLEX det = M00*M11 - M01*M10;
    if ( det == RZERO ){
        // fprintf(stderr, "%.5e+%.5ej %.5e+%.5ej \n", creal(M[0][0]), cimag(M[0][0]), creal(M[0][1]), cimag(M[0][1]));
        // fprintf(stderr, "%.5e+%.5ej %.5e+%.5ej \n", creal(M[1][0]), cimag(M[1][0]), creal(M[1][1]), cimag(M[1][1]));
        // fprintf(stderr, "matrix2x2 det=0.0, set matrix inv = 0.0.\n");
        // det = RZERO;
        *stats = INVERSE_FAILURE;
        return;
    }

    invM[0][0] = M11 / det;
    invM[0][1] = - M01 / det;
    invM[1][0] = - M10 / det;
    invM[1][1] = M00 / det;
    *stats = INVERSE_SUCCESS;
}

/**
 * 计算2x2复矩阵的和  
 * 
 * @param[in]      M1     矩阵1
 * @param[in]      M2     矩阵2
 * @param[out]     M      和矩阵
 */ 
inline GCC_ALWAYS_INLINE void cmat2x2_add(const MYCOMPLEX M1[2][2], const MYCOMPLEX M2[2][2], MYCOMPLEX M[2][2]){
    M[0][0] = M1[0][0] + M2[0][0];
    M[0][1] = M1[0][1] + M2[0][1];
    M[1][0] = M1[1][0] + M2[1][0];
    M[1][1] = M1[1][1] + M2[1][1];
}

/**
 * 计算2x2复矩阵的差  
 * 
 * @param[in]      M1     矩阵1
 * @param[in]      M2     矩阵2
 * @param[out]     M      差矩阵, M1 - M2
 */ 
inline GCC_ALWAYS_INLINE void cmat2x2_sub(const MYCOMPLEX M1[2][2], const MYCOMPLEX M2[2][2], MYCOMPLEX M[2][2]){
    M[0][0] = M1[0][0] - M2[0][0];
    M[0][1] = M1[0][1] - M2[0][1];
    M[1][0] = M1[1][0] - M2[1][0];
    M[1][1] = M1[1][1] - M2[1][1];
}

/**
 * 计算单位阵与2x2复矩阵的差  
 * 
 * @param[in,out]     M       差矩阵 I-M2
 */ 
inline GCC_ALWAYS_INLINE void cmat2x2_one_sub(MYCOMPLEX M[2][2]){
    M[0][0] = RONE - M[0][0];
    M[0][1] = - M[0][1];
    M[1][0] = - M[1][0];
    M[1][1] = RONE - M[1][1];
}

/**
 * 计算2x2复矩阵的积(矩阵相乘)  
 * 
 * @param[in]      M1     矩阵1
 * @param[in]      M2     矩阵2
 * @param[out]     M      积矩阵, M1 * M2
 */ 
inline GCC_ALWAYS_INLINE void cmat2x2_mul(const MYCOMPLEX M1[2][2], const MYCOMPLEX M2[2][2], MYCOMPLEX M[2][2]){
    MYCOMPLEX M011, M012, M021, M022;
    MYCOMPLEX M111, M112, M121, M122;
    M011 = M1[0][0]; M012 = M1[0][1]; 
    M021 = M1[1][0]; M022 = M1[1][1]; 
    M111 = M2[0][0]; M112 = M2[0][1]; 
    M121 = M2[1][0]; M122 = M2[1][1]; 
    M[0][0] = M011 * M111 + M012 * M121;
    M[0][1] = M011 * M112 + M012 * M122;
    M[1][0] = M021 * M111 + M022 * M121;
    M[1][1] = M021 * M112 + M022 * M122;
}

/**
 * 计算2x2复矩阵和常量的积
 * 
 * @param[in]      M1     矩阵1
 * @param[in]      k      常数
 * @param[out]     M      积矩阵, k * M2
 */
inline GCC_ALWAYS_INLINE void cmat2x2_k(const MYCOMPLEX M1[2][2], MYCOMPLEX k0, MYCOMPLEX M[2][2]){
    M[0][0] = M1[0][0] * k0;
    M[0][1] = M1[0][1] * k0;
    M[1][0] = M1[1][0] * k0;
    M[1][1] = M1[1][1] * k0;
}

/**
 * 计算2x2复矩阵和2x1的复向量的积
 * 
 * @param[in]      M1     矩阵1
 * @param[in]      M2     向量
 * @param[out]     M      积矩阵, M1 * M2
 */
inline GCC_ALWAYS_INLINE void cmat2x1_mul(const MYCOMPLEX M1[2][2], const MYCOMPLEX M2[2], MYCOMPLEX M[2]){
    MYCOMPLEX M00, M10;
    M00 = M1[0][0]*M2[0] + M1[0][1]*M2[1];
    M10 = M1[1][0]*M2[0] + M1[1][1]*M2[1];
    M[0] = M00;
    M[1] = M10;
}

/** 
 * 2x2复矩阵赋值 
 * 
 * @param[in]      M1     源矩阵
 * @param[out]     M2     目标矩阵
 */
inline GCC_ALWAYS_INLINE void cmat2x2_assign(const MYCOMPLEX M1[2][2], MYCOMPLEX M2[2][2]){
    M2[0][0] = M1[0][0];
    M2[0][1] = M1[0][1];
    M2[1][0] = M1[1][0];
    M2[1][1] = M1[1][1];
}

/** 
 * 计算nxn复矩阵的积(小矩阵)(最暴力的方式)
 * 
 * @param[in]     m1            M1矩阵行数
 * @param[in]     n1            M1矩阵列数
 * @param[in]     p1            M2矩阵列数
 * @param[in]     M1            M1矩阵 
 * @param[in]     M2            M2矩阵 
 * @param[out]    M             积矩阵 M1 * M2
 */
inline GCC_ALWAYS_INLINE void cmatmxn_mul(MYINT m1, MYINT n1, MYINT p1, const MYCOMPLEX M1[m1][n1], const MYCOMPLEX M2[n1][p1], MYCOMPLEX M[m1][p1]){
    MYINT m, n, k;
    MYCOMPLEX M0[m1][p1];
    for(m=0; m<m1; ++m){
        for(n=0; n<p1; ++n){
            M0[m][n] = 0.0;
            for(k=0; k<n1; ++k){
                M0[m][n] += M1[m][k] * M2[k][n];
            }
        }
    }

    // memcpy(M, M0, sizeof(MYCOMPLEX)*m1*p1);
    for(m=0; m<m1; ++m){
        for(n=0; n<p1; ++n){
            M[m][n] = M0[m][n];
        }
    }
}

/** 
 * 从M1大矩阵中划分Q子矩阵
 * 
 * @param[in]     m1           M1矩阵行数
 * @param[in]     n1           M1矩阵列数
 * @param[in]     M1           M1矩阵 
 * @param[in]     im           子矩阵起始行索引
 * @param[in]     in           子矩阵起始列索引
 * @param[in]     lm           子矩阵行数
 * @param[in]     ln           子矩阵列数
 * @param[out]    Q            子矩阵
 */
inline GCC_ALWAYS_INLINE void cmatmxn_block(MYINT m1, MYINT n1, const MYCOMPLEX M[m1][n1], MYINT im, MYINT in, MYINT lm, MYINT ln, MYCOMPLEX Q[lm][ln]){
    for(MYINT m=0; m<lm; ++m){
        for(MYINT n=0; n<ln; ++n){
            Q[m][n] = M[im+m][in+n];
        }
    }
}

/**
 * 打印矩阵 
 * 
 * @param[in]     m1          M1矩阵行数
 * @param[in]     n1          M1矩阵列数
 * @param[in]     M1          M1矩阵 
 * 
 */
inline GCC_ALWAYS_INLINE void cmatmxn_print(MYINT m1, MYINT n1, const MYCOMPLEX M1[m1][n1]){
    for(MYINT i=0; i<m1; ++i){
        for(MYINT j=0; j<n1; ++j){
            fprintf(stderr, " %15.5e + J%-15.5e ", creal(M1[i][j]), cimag(M1[i][j]));
        }
        fprintf(stderr, "\n");
    }
}