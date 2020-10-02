//
//  mynla.c
//  asd_sphere_20200508
//
//  Created by Mingyang Tan on 5/10/20.
//  Copyright © 2020 Mingyang Tan. All rights reserved.
//

#include "mynla.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 Find the [m-by-m] part of B that is the inverse of the [m-by-m] part of A, using Cholesky decomposition
 Input:
 m                   dimension of the matrix
 A                   matrix to be inverted
 lda                 leading order of A
 ldb                 leading order of B
 Output:
 B                   inverse of A
 */
void cholesky_inverse(const int m,
                      const double * A, const int lda,
                      double * B, const int ldb)
{
    int i, j, k;
    
    double s;
    
    double * l;
    double * linv;
    double * ident;
    
    l = calloc(m * m, sizeof(double));
    linv = calloc(m * m, sizeof(double));
    ident = calloc(m * m, sizeof(double));
    
    for (i = 0; i < m; i++) {
        ident[i * m + i] = 1.0;
    }
    
    /* Find Cholesky decomposition of A */
    for (i = 0; i < m; i++) {
        
        for (j = 0; j < m; j++) {
            
            if (i == j) {
                s = 0.0;
                for (k = 0; k < j; k++) {
                    s += l[j * m + k] * l[j * m + k];
                } // for k
                l[i * m + j] = sqrt(A[i * lda + j] - s);
                
            } else {
                s = 0.0;
                for (k = 0; k < j; k++) {
                    s += l[i * m + k] * l[j * m + k];
                }
                l[i * m + j] = 1.0 / l[j * m + j] * (A[i * lda + j] - s);
                
            } // if ... else ...
            
        } // for j
        
    } // for i
    
    /* Inverse of l */
    for (i = 0; i < m; i++) {
        
        for (j = i; j < m; j++) {
            s = 0.0;
            for (k = 0; k < j; k++) {
                s += l[j * m + k] * linv[k * m + i];
            } // for k
            linv[j * m + i] = 1.0 / l[j * m + j] * (ident[j * m + i] - s);
            
        } // for j
        
    } // for i
    
    /* B = (l^{T})^{-1} . l^{-1} */
    for (i = 0; i < m; i++) {
        
        for (j = 0; j < m; j++) {
            
            s = 0.0;
            for (k = 0; k < m; k++) {
                s += linv[k * m + j] * linv[k * m + i];
            }
            
            B[i * ldb + j] = s;
            
        } // for j
        
    } // for i
    
    free(l);
    free(linv);
    free(ident);
}

/*
 LU decomposition of the [m-by-m] part of A to obtain a [m-by-m] lower triangular of L and upper triangular of U
 Input:
 m                   dimension of the matrix
 A                   input matrix
 lda                 leading order of A
 ldl                 leading order of L
 ldu                 leading order of U
 Output:
 L                   lower-triangular matrix
 U                   upper-triangular matrix
 */
void lu_decomp(const int m,
               const double * A, const int lda,
               double * L, const int ldl,
               double * U, const int ldu)
{
    int i, j, k;
    
    double s;
    
    for (i = 0; i < m; i++) {
        // Upper triangular
        for (j = i; j < m; j++) {
            s = 0.0;
            for (k = 0; k < i; k++) {
                s += U[k * ldu + j] * L[i * ldl + k];
            } // for k
            U[i * ldu + j] = A[i * lda + j] - s;
        } // for j
        
        // Lower triangular
        for (j = i; j < m; j++) {
            if (j == i) {
                L[i * ldl + i] = 1.0;
            } else {
                s = 0.0;
                for (k = 0; k < i; k++) {
                    s += U[k * ldu + i] * L[j * ldl + k];
                }
                L[j * ldl + i] = (A[j * lda + i] - s) / U[i * ldl + i];
                
            } // if ... else ...
            
        } // for j
        
    } // for i
    
}
