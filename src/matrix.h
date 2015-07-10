// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   matrix.h
 * Author: fabien
 *
 * Created on July 9, 2010, 6:24 PM
 */

#ifndef _MATRIX_H
#define    _MATRIX_H

#ifdef    __cplusplus
extern "C" {
#endif

struct matrix_ {
    int     r;   // number of rows
    int     c;   // number of columns
    double  **e; // elements
};

#define MAT(x,i,j) x->e[i][j]

typedef struct matrix_ * Matrix;

long int MatrixAllocated();
Matrix  *MatrixArrayNew(int k, int n, int m); // array of k matrices dimension (n,m)
void    MatrixArrayDelete(Matrix *a, int k);
void    MatrixAxpy(Matrix y, double a, Matrix x);
Matrix  MatrixCopy(Matrix m);
void    MatrixDelete(Matrix m);
Matrix  MatrixEye(int n);
void    MatrixGadd(Matrix c, Matrix a, Matrix b, double alpha, double beta);
void    MatrixGemm(Matrix c, Matrix a, Matrix b, double alpha, double beta);
void    MatrixGemv(Vector y, Matrix A, Vector x, double alpha, double beta);
void    MatrixGemTv(Vector y, Matrix A, Vector x, double alpha, double beta);
double  MatrixGet(Matrix m, int i, int j);
void    MatrixGetSlice(Matrix A, int r0, int r1, int c0, int c1, Matrix B);
Matrix  MatrixNew(int r, int c);
double  MatrixMax(Matrix m);
double  MatrixMin(Matrix m);
void    MatrixMultTVec(Vector a, Matrix m, Vector c);
void    MatrixMultVec(Vector a, Matrix m, Vector c);
double  MatrixNorm(Matrix m);
double  MatrixNormInf(Matrix m);
double  MatrixNormOne(Matrix m);
void    MatrixPrint(Matrix m);
void    MatrixSet(Matrix m, int i, int j, double x);
void    MatrixSetAllTo(Matrix m, double x);
void    MatrixSetEqual(Matrix a, Matrix b);
void    MatrixScale(Matrix m, double s);
int     MatrixSolveLower(Matrix L, Vector b, Vector y);
int     MatrixSolveLowerT(Matrix L, Vector b, Vector y);
int     MatrixSolveUpper(Matrix U, Vector y, Vector x);
int     MatrixSolveUpperT(Matrix U, Vector y, Vector x);
void    MatrixTranspose(Matrix A, Matrix B);

#ifdef    __cplusplus
}
#endif

#endif    /* _MATRIX_H */
