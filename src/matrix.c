// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/*
 *  matrix.c
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "runtime.h"
#include "vector.h"
#include "matrix.h"

static long int _matrix_allocated_ = 0;

// y = a*x + y
void MatrixAxpy(Matrix y, double a, Matrix x) {
    int i, j;
    double *yi, *xi;

    if (BOUNDS_CHECK) {
        if ((y->r != x->r) || (y->c != x->c)) {
            RuntimeError("MatrixAxpy: invalid dimensions");
        }
    }
    for (i = 0; i < y->r; i++) {
        yi = y->e[i];
        xi = x->e[i];
        for (j = 0; j < y->c; j++) {
            yi[j] += a * xi[j];
        }
    }
}

// a = copy(b)
Matrix MatrixCopy(Matrix a) {
    Matrix b = MatrixNew(a->r, a->c);
    MatrixSetEqual(b, a);
    return b;
}

long int MatrixAllocated() {
    return (_matrix_allocated_);
}

void MatrixDelete(Matrix m) {
    int i;
    if (m == NULL) {
        return;
    }
    for (i = 0; i < m->r; i++) {
        free(m->e[i]);
    }
    free(m->e);
    free(m);
    _matrix_allocated_--;
}

// return a n-by-n identity matrix
Matrix MatrixEye(int n) {
    Matrix I;
    int i;
    if (n < 1) {
        RuntimeError("MatrixEye: n < 1");
    }
    I = MatrixNew(n, n);
    for (i = 0; i < n; i++) {
        I->e[i][i] = 1.0;
    }
    return I;
}

// c = alpha * a + beta * b
void MatrixGadd(Matrix c, Matrix a, Matrix b, double alpha, double beta) {
    int n = c->r;
    int m = c->c;
    int i, j;
    double *ci, *ai, *bi;
    if (BOUNDS_CHECK)
        if ((n != a->r) || (n != b->r) || (m != a->c) || (m != b->c)) {
            RuntimeError("MatrixGadd: invalid matrix dimensions");
        }
    for (i = 0; i < n; i++) { // do in parallel
        ci = c->e[i];
        ai = a->e[i];
        bi = b->e[i];
        for (j = 0; j < m; j++) {
            ci[j] = alpha * ai[j] + beta * bi[j];
        }
    }
}

// x = m->e[i][j]
double MatrixGet(Matrix m, int i, int j) {
    if (BOUNDS_CHECK) {
        if ((i < 0) || (i > m->r - 1)) {
            RuntimeError("MatrixGet: row index out of bounds");
        }
        if ((j < 0) || (j > m->c - 1)) {
            RuntimeError("MatrixGet: column index out of bounds");
        }
    }
    return m->e[i][j];
}

// B = A[r0:r1, c0:c1]
void MatrixGetSlice(Matrix A, int r0, int r1, int c0, int c1,
        Matrix B) {
    int n = A->r;
    int m = A->c;
    int nB = B->r;
    int mB = B->c;
    int i, j;
    if (BOUNDS_CHECK)
        if ((r0 > r1) || (c0 > c1) || (r0 < 0) || (c0 < 0) ||
                (r1 > n - 1) || (c1 > m - 1) || (nB != (r1 - r0 + 1)) || (mB != (c1 - c0 + 1))) {
            RuntimeError("MatrixGetSlice: invalid slice dimensions");
        }
    for (i = r0; i <= r1; i++) {
        for (j = c0; j <= c1; j++) {
            B->e[i - r0][j - c0] = A->e[i][j];
        }
    }
}

// c = alpha * a * b + beta * c
void MatrixGemm(Matrix c, Matrix a, Matrix b, double alpha, double beta) {
    int i, j, k;
    double s;
    double *ai, *ci;
    if (BOUNDS_CHECK)
        if ((c->r != a->r) || (c->c != b->c) || (a->c != b->r)) {
            RuntimeError("MatrixGemm: incompatible matrix dimensions");
        }

    for (i = 0; i < a->r; i++) {
        ai = a->e[i];
        ci = c->e[i];
        for (j = 0; j < b->c; j++) {
            s = 0.0;
            for (k = 0; k < a->c; k++) {
                s += ai[k] * b->e[k][j];
            }
            // TODO: test alpha = 1, beta = 0
            ci[j] = alpha * s + beta * ci[j];
        }
    }
}

// y = alpha * A * x + beta * y
void MatrixGemv(Vector y, Matrix A, Vector x, double alpha, double beta) {
    int i, j;
    double s;
    double *Ai;
    if (BOUNDS_CHECK)
        if ((y->r != A->r) || (A->c != x->r)) {
            RuntimeError("MatrixGemv: incompatible matrix dimensions");
        }
    for (i = 0; i < A->r; i++) {
        s = 0.0;
        Ai = A->e[i];
        for (j = 0; j < A->c; j++) {
            s += Ai[j] * x->e[j];
        }
        y->e[i] = alpha * s + beta * y->e[i];
    }
}

// y = alpha * A^T * x + beta * y
void MatrixGemTv(Vector y, Matrix A, Vector x, double alpha, double beta) {
    int i, j;
    double s;
    if (BOUNDS_CHECK)
        if ((y->r != A->c) || (A->r != x->r)) {
            RuntimeError("MatrixGemTv: incompatible matrix dimensions");
        }
    for (i = 0; i < A->c; i++) {
        s = 0.0;
        for (j = 0; j < A->r; j++) {
            s += A->e[j][i] * x->e[j];
        }
        y->e[i] = alpha * s + beta * y->e[i];
    }
}

// max(m->e[i][j])
double MatrixMax(Matrix m) {
    int i, j;
    double *mi;
    double x = m->e[0][0];
    for (i = 0; i < m->r; i++) {
        mi = m->e[i];
        for (j = 0; j < m->c; j++) {
            if (mi[j] > x) {
                x = mi[j];
            }
        }
    }
    return x;
}

// min(m->e[i][j])
double MatrixMin(Matrix m) {
    int i, j;
    double *mi;
    double x = m->e[0][0];
    for (i = 0; i < m->r; i++) {
        mi = m->e[i];
        for (j = 0; j < m->c; j++) {
            if (mi[j] < x) {
                x = mi[j];
            }
        }
    }
    return x;
}

// a = M*c
void MatrixMultVec(Vector a, Matrix m, Vector c) {
    int i, j;
    double x;
    double *mi;
    if (BOUNDS_CHECK)
        if ((a->r != m->r) || (m->c != c->r)) {
            RuntimeError("MatrixMultVec: invalid dimensions");
        }
    for (i = 0; i < a->r; i++) {
        x = 0.0;
        mi = m->e[i];
        for (j = 0; j < m->c; j++) {
            x += mi[j] * c->e[j];
        }
        a->e[i] = x;
    }
}

// a = M^T*c
void MatrixMultTVec(Vector a, Matrix m, Vector c) {
    int i, j;
    double x;
    if (BOUNDS_CHECK)
        if ((a->r != m->c) || (m->r != c->r)) {
            RuntimeError("MatrixMultTVec: invalid dimensions");
        }
    for (i = 0; i < a->r; i++) {
        x = 0.0;
        for (j = 0; j < c->r; j++) {
            x += m->e[j][i] * c->e[j];
        }
        a->e[i] = x;
    }
}

// create a new r-by-c matrix
Matrix MatrixNew(int r, int c) {
    int i;
    Matrix m;
    if (r < 1) {
        RuntimeError("Matrix_new: row dimension < 1");
    }
    if (c < 1) {
        RuntimeError("Matrix_new: column dimension < 1");
    }
    m = (Matrix) malloc(sizeof (struct matrix_));
    m->r = r;
    m->c = c;
    m->e = (double **) malloc(r * sizeof (double *));
    for (i = 0; i < r; i++) {
        m->e[i] = (double *) calloc(c, sizeof (double));
    }
    _matrix_allocated_++;
    return m;
}

// Construct and array of k matricies with dimension (n,m)
Matrix *MatrixArrayNew(int k, int n, int m) {
    int i;
    Matrix *a;
    if ((k < 1) || (n < 1) || (m < 1)) {
        RuntimeError("MatrixArray: invalid dimensions");
    }
    a = (Matrix *) malloc(k * sizeof (Matrix));
    for (i = 0; i < k; i++) {
        a[i] = MatrixNew(n, m);
    }
    return a;
}

// delete the array of matrices
void MatrixArrayDelete(Matrix *a, int k) {
    int i;
    if (a == NULL) {
        return;
    }
    for (i = 0; i < k; i++) {
        MatrixDelete(a[i]);
    }
    free(a);
}

// sqrt(sum(m->e[i][j]^2))
double MatrixNorm(Matrix m) {
    int i, j;
    double *mi;
    double x = 0.0;
    for (i = 0; i < m->r; i++) {
        mi = m->e[i];
        for (j = 0; j < m->c; j++) {
            x += mi[j] * mi[j];
        }
    }
    return sqrt(x);
}

// max_i(sum_j|m->e[i][j]|)
double MatrixNormInf(Matrix m) {
    int i, j;
    double rs;
    double *mi;
    double x = 0.0;
    for (i = 0; i < m->r; i++) {
        rs = 0.0;
        mi = m->e[i];
        for (j = 0; j < m->c; j++) {
            rs += fabs(mi[j]);
        }
        if (rs > x) {
            x = rs;
        }
    }
    return x;
}

// max_j(sum_i|m->e[i][j]|)
double MatrixNormOne(Matrix m) {
    int i, j;
    double cs;
    double x = 0.0;
    for (j = 0; j < m->c; j++) {
        cs = 0.0;
        for (i = 0; i < m->r; i++) {
            cs += fabs(m->e[i][j]);
        }
        if (cs > x) {
            x = cs;
        }
    }
    return x;
}

// print a matrix to standard output
void MatrixPrint(Matrix m) {
    int i, j;
    if (m == NULL) {
        printf("NULL\n");
        return;
    }
    for (i = 0; i < m->r; i++) {
        for (j = 0; j < m->c; j++) {
            printf("%g ", m->e[i][j]);
        }
        printf("\n");
    }
}

// m = scalar*m
void MatrixScale(Matrix m, double s) {
    int i, j;
    double *mi;
    for (i = 0; i < m->r; i++) {
        mi = m->e[i];
        for (j = 0; j < m->c; j++) {
            mi[j] *= s;
        }
    }
}

// m->e[i][j] = x
void MatrixSet(Matrix m, int i, int j, double x) {
    if (BOUNDS_CHECK) {
        if ((i < 0) || (i > m->r - 1)) {
            RuntimeError("MatrixSet: row index out of bounds");
        }
        if ((j < 0) || (j > m->c - 1)) {
            RuntimeError("MatrixSet: column index out of bounds");
        }
    }
    m->e[i][j] = x;
}

// m->e[i][j] = x, for all i,j
void MatrixSetAllTo(Matrix m, double x) {
    int i, j;
    double *mi;
    for (i = 0; i < m->r; i++) {
        mi = m->e[i];
        for (j = 0; j < m->c; j++) {
            mi[j] = x;
        }
    }
}

// a = b
void MatrixSetEqual(Matrix a, Matrix b) {
    int i, j;
    double *ai, *bi;
    if (BOUNDS_CHECK)
        if ((a->r != b->r) || (a->c != b->c)) {
            RuntimeError("MatrixSetEqual: invalid dimensions");
        }
    for (i = 0; i < a->r; i++) {
        ai = a->e[i];
        bi = b->e[i];
        for (j = 0; j < a->c; j++) {
            ai[j] = bi[j];
        }
    }
}

// Solve a lower triangular system, i.e., L*y = b (forward substitution)
// L nonsingular lower triangular n-by-n Matrix
// b righthand side: n-Vector
// y solution: n-Vetor
// return 0 if successful, != 0 otherwise
int MatrixSolveLower(Matrix L, Vector b, Vector y) {
    int i, k;
    int n = L->r;
    double s;
    double *Li;

    if (BOUNDS_CHECK)
        if ((n != L->c) || (n != b->r) || (n != y->r)) {
            RuntimeWarning("SolveLower: incompatible dimensions for L, x, y");
            return 1;
        }

    for (i = 0; i < n; i++) {
        Li = L->e[i];
        if (fabs(Li[i]) < DBL_EPSILON) {
            //RuntimeWarning("SolveLower: |L[i][i]| < DBL_EPSILON");
            return 1;
        }
        s = 0.0;
        for (k = 0; k < i; k++) {
            s += Li[k] * y->e[k];
        }
        y->e[i] = (b->e[i] - s) / Li[i];
    }
    return 0;
}

// Solve a upper triangular system, (L^T)*y = b (backward substitution)
// L nonsingular lower triangular n-by-n matrix
// b righthand side: n-Vector
// y solution: n-Vector
// return 0 if successful, != 0 otherwise
int MatrixSolveLowerT(Matrix L, Vector b, Vector y) {
    int i, k;
    double s;
    int n = L->r;
    double *Li;
    if (BOUNDS_CHECK)
        if ((n != L->c) || (n != b->r) || (n != y->r)) {
            RuntimeWarning("SolveLowerT: incompatible dimensions for L, x, y");
            return 1;
        }
    for (i = n - 1; i >= 0; i--) {
        Li = L->e[i];
        if (fabs(Li[i]) < DBL_EPSILON) {
            //RuntimeWarning("SolveLowerT: |L[i][i]| < DBL_EPSILON");
            return 1;
        }
        s = 0.0;
        for (k = i + 1; k < n; k++) {
            s += L->e[k][i] * y->e[k];
        }
        y->e[i] = (b->e[i] - s) / Li[i];
    }
    return 0;
}

// Solve a upper triangular system, i.e., U*x = y (backward substitution)
// U nonsingular upper triangular n-by-n matrix
// y righthand side: n-Vector
// x solution: n-Vector
// return 0 if successful, != 0 otherwise
int MatrixSolveUpper(Matrix U, Vector y, Vector x) {
    int i, k;
    int n = U->r;
    double s;
    double *Ui;
    if (BOUNDS_CHECK)
        if ((n != U->c) || (n != y->r) || (n != x->r)) {
            RuntimeWarning("SolveUpper: incompatible dimensions for U, x, y");
            return 1;
        }
    for (i = n - 1; i >= 0; i--) {
        Ui = U->e[i];
        if (fabs(Ui[i]) < DBL_EPSILON) {
            //RuntimeWarning("SolveUpperT: |U[i][i]| < DBL_EPSILON");
            return 1;
        }
        s = 0.0;
        for (k = i + 1; k < n; k++) {
            s += Ui[k] * x->e[k];
        }
        x->e[i] = (y->e[i] - s) / Ui[i];
    }
    return 0;
}

// Solve U^T*x = y, U n-by-n upper triangular
// (same as Solve L*x = y)
int MatrixSolveUpperT(Matrix U, Vector y, Vector x) {
    double s;
    int n = U->r;
    int i, k;
    double *Ui;
    //n, m := U.Size()
    if (BOUNDS_CHECK)
        if ((n != U->c) || (n != y->r) || (n != x->r)) {
            RuntimeWarning("SolveUpperT: incompatible dimensions for U, x, y");
            return 2;
        }
    for (i = 0; i < n; i++) {
        Ui = U->e[i];
        if (fabs(Ui[i]) < DBL_EPSILON) {
            //RuntimeWarning("SolveUpperT: |U[i][i]| < DBL_EPSILON");
            return 1;
        }
        s = 0.0;
        for (k = 0; k < i; k++) {
            s += U->e[k][i] * x->e[k];
        }
        x->e[i] = (y->e[i] - s) / Ui[i];
    }

    return 0;
}

// A = B^T
void MatrixTranspose(Matrix A, Matrix B) {
    int i, j;
    double *Ai;
    if (BOUNDS_CHECK)
        if ((A->r != B->c) || (B->r != A->c)) {
            RuntimeError("MatrixTranspose: invalid dimensions");
        }
    for (i = 0; i < A->r; i++) {
        Ai = A->e[i];
        for (j = 0; j < A->c; j++) {
            Ai[j] = B->e[j][i];
        }
    }
}
