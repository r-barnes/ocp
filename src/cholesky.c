// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "runtime.h"
#include "vector.h"
#include "matrix.h"
#include "cholesky.h"

// Factor a symmetric positive definite matrix.
// A symmetic positive definiet matrix, n-by-n Matrix
// L lower triangular matrix such that A = L*L^T, n-by-n Matrix
// return 0 normal completion; 1 singular matrix/not positive definite.
int CholeskyFactor(Matrix A, Matrix L) {
    int i, j, k;
    double s;
    double *Li, *Lj, *Ai;

    int n = A->r;

    if ((n != A->c) || (n != L->r) || (n != L->c)) {
        RuntimeWarning("CholeskyFactor: invalid dimensions");
        return 1;
    }

    MatrixSetAllTo(L, 0.0);

    for (i = 0; i < n; i++) {
        Li = L->e[i];
        Ai = A->e[i];
        for (j = 0; j < i; j++) {
            s = 0.0;
            Lj = L->e[j];
            for (k = 0; k < j; k++) {
                s += Li[k] * Lj[k];
            }
            Li[j] = (Ai[j] - s) / Lj[j];
        }
        s = 0;
        for (k = 0; k < i; k++) {
            s += Li[k] * Li[k];
        }
        s = Ai[i] - s;
        if (s > DBL_EPSILON) {
            Li[i] = sqrt(s);
        } else {
            return 1;
        }
    }
    return 0;
}

// Solve a symmetric positive definite linear system, A*x = b.
// A symmetic positive definiet matrix, n-by-n Matrix
// b righthand side of the equations, n-Vector
// x solution, n-Vector
// return 0 normal completion; 1 singular matrix/not positive definite.
int CholeskySolveLinearSystem(Matrix A, Vector b, Vector x) {
    int error;
    Matrix L;
    Vector y;
    int n = A->r;
    if ((n != A->c) || (n != b->r) || (n != x->r)) {
        RuntimeWarning("CholeskySolveLinearSystem: invalid dimensions");
        return 1;
    }
    L = MatrixNew(n, n);
    error = CholeskyFactor(A, L);
    if (error != 0) {
        RuntimeWarning("CholeskySolveLinearSystem: matrix not positive definite");
        MatrixDelete(L);
        return 1;
    }
    y = VectorNew(n);
    error = MatrixSolveLower(L, b, y);
    if (error != 0) {
        RuntimeWarning("Cholesky error: singular matrix");
        MatrixDelete(L);
        VectorDelete(y);
        return 1;
    }
    error = MatrixSolveLowerT(L, y, x);
    if (error != 0) {
        RuntimeWarning("Cholesky error: singular matrix");
        MatrixDelete(L);
        VectorDelete(y);
        return 1;
    }
    MatrixDelete(L);
    VectorDelete(y);
    return 0;
}
