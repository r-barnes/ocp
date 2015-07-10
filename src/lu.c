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
#include "lu.h"

int MinInt(int i, int j) {
    int r = (i < j ? i : j);
    return r;
}

// Solve L*U*x = P*y, P[i][p[i]] = 1, otherwise P[i][j] = 0
int LUSolve(Matrix L, Matrix U, int *p, Vector y, Vector x) {
    if (BOUNDS_CHECK) {
        if ((L->r != L->c) || (L->r != U->r) || (L->r != U->c) ||
                (L->r != y->r) || (L->r != x->r)) {
            RuntimeWarning("LUSolve: incompatible matrix dimensions");
            return 1;
        }
    }
    
    int n = L->r;
    Vector c, v;
    int i, err;
    c = VectorNew(n);
    v = VectorNew(n);

    for (i = 0; i < n; i++) {
        c->e[i] = y->e[p[i]];
    }

    err = MatrixSolveLower(L, c, v);
    if (err != 0) {
        VectorDelete(c);
        VectorDelete(v);
        return 1;
    }
    
    err = MatrixSolveUpper(U, v, x);
    if (err != 0) {
        VectorDelete(c);
        VectorDelete(v);
        return 1;
    }

    VectorDelete(c);
    VectorDelete(v);
    return 0;
}

// Solve L*U*X = P*Y, P[i][p[i]] = 1, otherwise P[i][j] = 0
int LUSolveM(Matrix L, Matrix U, int *p, Matrix Y, Matrix X) {
    if (BOUNDS_CHECK) {
        if ((L->r != L->c) || (L->r != U->r) || (L->r != U->c) ||
                (L->r != Y->r) || (L->r != X->r) ||
                (Y->c != X->c)) {
            RuntimeWarning("LUSolve: incompatible matrix dimensions");
            return 1;
        }
    }
    int i, j, n, m, err = 0;
    Vector x, y;
    n = L->r;
    m = Y->c;
    x = VectorNew(n);
    y = VectorNew(n);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            y->e[j] = Y->e[j][i];
        }
        err = LUSolve(L, U, p, y, x);
        if (err != 0) {
            break;
        }
        for (j = 0; j < n; j++) {
            X->e[j][i] = x->e[j];
        }
    }
    VectorDelete(x);
    VectorDelete(y);
    return err;
}

// Returns the LU factors of P * M, i.e., P * M = L * U
// -: err == 0, successful completion
// -: err == 1 if M is singular, i.e., the absolute value
//   of some diagonal element of U is <= DoubleEpsilon
// -: err == 2 if the dimensions of M, L, U, p are incompatible
// -: M.Size() = (n, n)
// -: L.Size() = (n, n)
// -: U.Size() = (n, m)
// -: len(p) = n
// -: P.Size() = (n, n)
// -: P[i][p[i]] = 1
int LUFactor(Matrix M, Matrix L, Matrix U, int *p) {
    int n = M->r;
    if (BOUNDS_CHECK) {
        if ((n != M->c) || (n != L->r) || (n != L->c) || (n != U->r) || (n != U->c)) {
            RuntimeWarning("LUFactor: incompatible matrix dimensions");
            return 1;
        }
    }
    int i, j, k;
    Matrix A = MatrixCopy(M);
    MatrixSetAllTo(U,0.0);
    MatrixSetAllTo(L,0.0);
    for (i = 0; i < n; i++) {
        p[i] = i;
    }
	for (j = 1; j <= n; j++) {
		for (i = 1; i <= j-1; i++) {
			double s = 0.0;
			for (k = 1; k <= i-1; k++) {
				s += L->e[i-1][k-1]*U->e[k-1][j-1];
			}
			U->e[i-1][j-1] = A->e[i-1][j-1] - s;
		}
		int r = j;
		double ur = 0.0;
		for (i = j; i <= n; i++) {
			double s = 0.0;
			for (k = 1; k <= j-1; k++) {
				s += L->e[i-1][k-1]*U->e[k-1][j-1];
			}
			U->e[i-1][j-1] = A->e[i-1][j-1] - s;
			if (fabs(U->e[i-1][j-1]) > ur) {
				ur = fabs(U->e[i-1][j-1]);
				r = i;
			}
		}
		if (ur < DBL_EPSILON) {
			MatrixDelete(A);
			return -1;
		}
		if (r != j) {
			int pj = p[j-1]; p[j-1] = p[r-1]; p[r-1] = pj;
			double *t;
			t = A->e[j-1]; A->e[j-1] = A->e[r-1]; A->e[r-1] = t;
			t = L->e[j-1]; L->e[j-1] = L->e[r-1]; L->e[r-1] = t;
			t = U->e[j-1]; U->e[j-1] = U->e[r-1]; U->e[r-1] = t;
		}
		for (i = j+1; i <= n; i++) {
			L->e[i-1][j-1] = U->e[i-1][j-1] / U->e[j-1][j-1];
		}
		L->e[j-1][j-1] = 1.0;
	}
	MatrixDelete(A);
	return 0;
}

/*
int LUFactor0(Matrix M, Matrix L, Matrix U, int *p) {
    int n = M->r;

    if (BOUNDS_CHECK) {
        if ((n != M->c) || (n != L->r) || (n != L->c) || (n != U->r) || (n != U->c)) {
            RuntimeWarning("LUFactor: incompatible matrix dimensions");
            return 1;
        }
    }
    
    int i, j, k;

    Matrix A = MatrixCopy(M);
    MatrixSetAllTo(U,0.0);
    MatrixSetAllTo(L,0.0);

    for (i = 0; i < n; i++) {
        p[i] = i;
    }
	
	Vector vv = VectorNew(n);
	
	for (i = 0; i < n; i++) {
		double amax = 0.0;
		for (j = 0; j < n; j++) {
			if (fabs(A->e[i][j]) > amax) {
				amax = fabs(A->e[i][j]);
			}
		}
		if (amax == 0.0) {
			VectorDelete(vv);
			MatrixDelete(A);
			return -1;
		}
		vv->e[i] = 1.0/amax;
	}
    
    for (j = 0; j < n; j++) {
    	for (i = 0; i <= j-1; i++) {
    		double s = A->e[i][j];
    		for (k = 0; k <= i-1; k++) {
    			s -= A->e[i][k]*A->e[k][j];
    		}
    		A->e[i][j] = s;
    	}
		double amax  = 0.0;
		int imax = j;
		for (i = j; i < n; i++) {
			double s = A->e[i][j];
			for (k = 0; k <= j-1; k++) {
				s -= A->e[i][k]*A->e[k][j];
			}
			A->e[i][j] = s;
			double dd = vv->e[i]*fabs(s);
			if (dd >= amax) {
				imax = i;
				amax = dd;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				double dd = A->e[imax][k];
				A->e[imax][k] = A->e[j][k];
				A->e[j][k] = dd;
			}
			vv->e[imax] = vv->e[j];
			int pj = p[j];
			p[j] = p[imax];
			p[imax] = pj;
		}
		
		if (fabs(A->e[j][j]) < DBL_EPSILON) {
			VectorDelete(vv);
			MatrixDelete(A);
			return -1;
		}
		
		if (j != n-1) {
			double dd = 1.0/A->e[j][j];
			for (i = j+1; i < n; i++) {
				A->e[i][j] *= dd;
			}
		}
    }
    
    for (i = 0; i < n; i++) {
    	for (j = 0; j < i; j++) {
    		L->e[i][j] = A->e[i][j];
    	}
    	L->e[i][i] = 1.0;
    	for (j = i; j < n; j++) {
    		U->e[i][j] = A->e[i][j];
    	}
    }
    MatrixDelete(A);
    VectorDelete(vv);
    return 0;
}


int LUFactor1(Matrix M, Matrix L, Matrix U, int *p) {
    int n = M->r;

    if (BOUNDS_CHECK) {
        if ((n != M->c) || (n != L->r) || (n != L->c) || (n != U->r) || (n != U->c)) {
            RuntimeWarning("LUFactor: incompatible matrix dimensions");
            return 1;
        }
    }
    
    int i, j, k;

    Matrix A = MatrixCopy(M);
    MatrixSetAllTo(U,0.0);
    MatrixSetAllTo(L,0.0);

    for (i = 0; i < n; i++) {
        p[i] = i;
    }
    // TODO: do this without l and u
    Vector l = VectorNew(n);
    Vector u = VectorNew(n);
    
    for (j = 0; j < n; j++) {
        double a = fabs(A->e[j][j]);
        int r = j;
        for (i = j+1; i < n; i++) {
            double b = fabs(A->e[i][j]);
            if (b > a) {
                a = b;
                r = i;
            }
        }
        if (r != j) { // swap rows r and j
            double *tmp;
            tmp = A->e[r];
            A->e[r] = A->e[j];
            A->e[j] = tmp;
            int t;
            t = p[r];
            p[r] = p[j];
            p[j] = t;
        }
        double lu = 0.0;
        for (i = 0; i < j; i++) {
            double su = 0;
            double sl = 0;
            for (k = 0; k < i; k++) {
                su = su + L->e[i][k] * u->e[k];
                sl = sl + U->e[k][i] * l->e[k];
            }
            u->e[i] = (A->e[i][j] - su);
            l->e[i] = (A->e[j][i] - sl) / U->e[i][i];
            lu = lu + l->e[i] * u->e[i];
            
            L->e[j][i] = l->e[i];
            U->e[i][j] = u->e[i];
        }
        L->e[j][j] = 1.0;
        U->e[j][j] = A->e[j][j] - lu; 
        if (fabs(U->e[j][j]) < DBL_EPSILON) {
            MatrixDelete(A);
            VectorDelete(l);
            VectorDelete(u);
            return 1;
        }
    }
    MatrixDelete(A);
    VectorDelete(l);
    VectorDelete(u);
    return 0;
}
*/

// Solve the linear system A*x = b
int LUSolveLinearSystem(Matrix A, Vector b, Vector x) {
	
	if ((A->r != A->c) || (A->r != b->r) || (A->r != x->r)) {
		RuntimeWarning("LUSolveLinearSystem(): error, invalid system dimensions");
		return -1;
	}
	int n = A->r;
	Matrix L = MatrixNew(n, n);
	Matrix U = MatrixNew(n, n);
	int *p = (int *) calloc(n, sizeof (int));
	
	int err = LUFactor(A, L, U, p);
	
	if (err != 0) {
		RuntimeWarning("LUSolveLinearSystem(): unable to factor matrix");
		MatrixDelete(L);
		MatrixDelete(U);
		free(p);
		return -2;
	}

	err = LUSolve(L, U, p, b, x);

	if (err != 0) {
		RuntimeWarning("LUSolveLinearSystem(): unable to solve");
		MatrixDelete(L);
		MatrixDelete(U);
		free(p);
		return -3;
	}	
	
	MatrixDelete(L);
	MatrixDelete(U);
	free(p);
	
	return 0;
}
