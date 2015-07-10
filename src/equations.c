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
#include "equations.h"

Equations EquationsNew() {
    Equations eqns;
    eqns = (Equations) malloc(sizeof (struct equations_));
    eqns->F = NULL;
    eqns->DF = NULL;
    // private
    eqns->x1 = NULL;
    eqns->F0 = NULL;
    eqns->F1 = NULL;
    return eqns;
}

void EquationsDelete(Equations eqns) {
    if (eqns->x1 != NULL) {
        VectorDelete(eqns->x1);
    }
    if (eqns->F0 != NULL) {
        VectorDelete(eqns->F0);
    }
    if (eqns->F1 != NULL) {
        VectorDelete(eqns->F1);
    }
    free(eqns);
}

// Forward difference Jacobian approximation
// Df: m-by-n Matrix: ~= df/dx
// x0: n-Vector
void EquationsDF(Equations eqns, Vector x0, Matrix Df) {
    int n = x0->r;
    int m = Df->r;
    int i, j;
    double delta;
    Vector x1, F0, F1;

    if (Df->c != n) {
        RuntimeWarning("EquationsJacobian: invalid dimensions.");
        return;
    }

    if (eqns->x1 == NULL) {
        eqns->x1 = VectorNew(n);
        eqns->F0 = VectorNew(m);
        eqns->F1 = VectorNew(m);
    }

    x1 = eqns->x1;
    F0 = eqns->F0;
    F1 = eqns->F1;

    VectorSetEqual(x1, x0);
    eqns->F(x0, F0);
    for (i = 0; i < n; i++) {
        delta = SQRT_EPS * fmax(1.0, fabs(x0->e[i]));
        x1->e[i] += delta;
        delta = x1->e[i] - x0->e[i];
        eqns->F(x1, F1);
        for (j = 0; j < m; j++) {
            Df->e[j][i] = (F1->e[j] - F0->e[j]) / delta;
        }
        x1->e[i] = x0->e[i];
    }
}
