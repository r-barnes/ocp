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
#include "qr.h"

void _householder(double **Q, double **R, double **B, int n, int m) {
    int i, j, k;
    Vector Beta = VectorNew(m);
    Vector W = VectorNew(n);
    Matrix V = MatrixNew(n, n);
    double *beta = Beta->e;

    // Q = I
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                Q[i][j] = 1.0;
            } else {
                Q[i][j] = 0.0;
            }
        }
    }
    
    // R = B
    for (i = 0; i < n; i++) {
        double *Ri = R[i];
        double *Bi = B[i];
        for (k = 0; k < m; k++) {
            Ri[k] = Bi[k];
        }
    }

    // find the Householder vectors for each column of R
    for (k = 0; k < m; k++) {
        double xi;
        double norm_x = 0.0;
        double *w = W->e;
        double *v = V->e[k];
        for (i = k; i < n; i++) {
            xi = R[i][k];
            norm_x += xi * xi;
            v[i] = xi;
        }
        norm_x = sqrt(norm_x);
        double x1 = v[k];
        double sign_x1 = sign(x1);
        double gamma = -1.0 * sign_x1 * norm_x;
        double v_dot_v = 2.0 * norm_x * (norm_x + sign_x1 * x1);

        if (v_dot_v < DBL_EPSILON) {
            for (i = k; i < m; i++) {
                R[i][i] = 0.0;
            }
            VectorDelete(Beta);
            VectorDelete(W);
            MatrixDelete(V);
            return;
        }
        v[k] -= gamma;
        beta[k] = -2.0 / v_dot_v;
      
        // w(k:m) = R(k:n, k:m)^T * v(k:n)
        for (i = k; i < m; i++) {
            double s = 0.0;
            for (j = k; j < n; j++) { // FIX: make this row-wise
                s += R[j][i] * v[j];
            }
            w[i] = s;
        }

        // R(k:n, k:m) += beta * v(k:n) * w(k:m)^T
        //a[k : n, k : m] = a[k : n, k : m] + beta * (v * wT)
        for (i = k; i < n; i++) {
            for (j = k; j < m; j++) {
                R[i][j] += beta[k] * v[i] * w[j];
            }
        }
    }

    for (k = m-1; k >= 0; k--) {
        double *v = V->e[k];
        double *u = W->e;
        //uT = v.transpose() * Q[k : n, k : n]
        for (i = k; i < n; i++) {
            double s = 0.0;
            for (j = k; j < n; j++) {
                s += Q[j][i] * v[j];
            }
            u[i] = s;
        }
        //Q[k : n, k : n] = Q[k : n, k : n] + Beta[k] * (v * uT)
        for (i = k; i < n; i++) {
            for (j = k; j < n; j++) {
                Q[i][j] += beta[k] * v[i] * u[j];
            }
        }
    }
    VectorDelete(Beta);
    VectorDelete(W);
    MatrixDelete(V);
}
