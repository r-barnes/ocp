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
#include "nlp.h"

static NLP thisNLP;

void NLPfhg(Vector F, Vector x);

// create a new NLP data structure
//
// public members
//
// n: number of unknown variables
// m: number of equality constraints
// p: number of inequality constraints
// fhg: double (*fhg)(Vector x, Vector h, Vector g)
//      input: x: n-Vector  (unknown variables)
//      output: h: m-Vector (equality constraints)
//      output: g: p-Vector (inequality constraints)
// Dfhg: void (*Dfhg)(Vector x, Vector Df, Matrix Dh, Matrix Dg): (optional)
//      input: x: n-Vector  (unknown variables)
//      output: Df: n-Vector (function gradient)
//      output: Dh: m-by-n Matrix (Jacobian equality constraints)
//      output: Dg: p-by-n Matrix (Jacobian inequality constraints)
NLP NLPNew() {
    NLP a;
    a = (NLP)malloc(sizeof(struct nlp_));
    a->n = 0;
    a->m = 0;
    a->p = 0;
    a->fhg = NULL;
    a->Dfhg = NULL;
    // private data for finite difference derivatives
    a->e = NULL;
    a->De = NULL;
    a->h = NULL;
    a->g = NULL;
    return a;
}

// delete the NLP data structure
void NLPDelete(NLP a) {
    if (a->e != NULL) {
        EquationsDelete(a->e);
    }
    if (a->De != NULL) {
        MatrixDelete(a->De);
    }
    if (a->h != NULL) {
        VectorDelete(a->h);
    }
    if (a->g != NULL) {
        VectorDelete(a->g);
    }
    free(a);
}

// compute the derivatives using a forward difference approximation
void NLPDfhg(NLP a, Vector x0, Vector Df, Matrix Dh, Matrix Dg) {
    int n, m, p, i, j;
    double *Dhi, *Dgi, *Dei;
    thisNLP = a;
    n = a->n;
    m = a->m;
    p = a->p;
    if (a->e == NULL) {
        a->e = EquationsNew();
        a->e->F = NLPfhg;
        a->e->DF = NULL;
        a->De = MatrixNew(1+m+p, n);
        if (m > 0) {
            a->h = VectorNew(m);
        }
        if (p > 0) {
            a->g = VectorNew(p);
        }
    }

    EquationsDF(a->e, x0, a->De);

    for (i = 0; i < n; i++) {
        Df->e[i] = a->De->e[0][i];
    }
    for (i = 0; i < m; i++) {
        Dhi = Dh->e[i];
        Dei = a->De->e[1+i];
        for (j = 0; j < n; j++) {
            Dhi[j] = Dei[j];
        }
    }
    for (i = 0; i < p; i++) {
        Dgi = Dg->e[i];
        Dei = a->De->e[1+m+i];
        for (j = 0; j < n; j++) {
            Dgi[j] = Dei[j];
        }
    }
}

// evaluate F(x) = [f(x); h(x); g(x)]
void NLPfhg(Vector x, Vector F) {
    int i;
    NLP a = thisNLP;
    F->e[0] = a->fhg(x, a->h, a->g);
    for (i = 0; i < a->m; i++) {
        F->e[1+i] = a->h->e[i];
    }
    for (i = 0; i < a->p; i++) {
        F->e[1+a->m+i] = a->g->e[i];
    }
}
