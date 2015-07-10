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
#include "qp.h"

// create a new QP data structure
QP QPNew() {
    QP a;
    a = (QP) malloc(sizeof (struct qp_));
    a->H = NULL;
    a->c = NULL;
    a->A1 = NULL;
    a->b1 = NULL;
    a->A2 = NULL;
    a->b2 = NULL;
    return a;
}

// Returns int *a = {error_code, n, m, p}
// If error_code != 0 the problem specification is invalid
void QPcheck(QP q, int *a) {
    int n = 0, m = 0, p = 0;
    Matrix H = q->H, A1 = q->A1, A2 = q->A2;
    Vector c = q->c, b1 = q->b1, b2 = q->b2;

    a[0] = 1;
    if ((H == NULL) || (c == NULL)) {
        RuntimeWarning("QP: H or c is NULL");
        return;
    }
    n = H->r;
    if (n != H->c) {
        RuntimeWarning("QP: H not square");
        return;
    }
    if (n != c->r) {
        RuntimeWarning("QP: invalid dimension for c");
        return;
    }
    if (A1 != NULL) {
        if (b1 == NULL) {
            RuntimeWarning("QP: b1 == NULL");
            return;
        }
        if (n != A1->c) {
            RuntimeWarning("QP: invalid dimension for A1");
            return;
        }
        m = A1->r;
        if (m != b1->r) {
            RuntimeWarning("QP: invalid dimension for b1");
            return;
        }
    }
    if (A2 != NULL) {
        if (b2 == NULL) {
            RuntimeWarning("QP: b2 == NULL");
            return;
        }
        if (n != A2->c) {
            RuntimeWarning("QP: invalid dimension for A2");
            return;
        }
        p = A2->r;
        if (p != b2->r) {
            RuntimeWarning("QP: invalid dimension for b2");
            return;
        }
    }
    a[0] = 0;
    a[1] = n;
    a[2] = m;
    a[3] = p;
}

void QPDelete(QP a) {
    if (a == NULL) {
        return;
    }
    free(a);
}

void QPPrint(QP a) {
    printf("H:\n");
    if (a->H != NULL) MatrixPrint(a->H);
    printf("c: ");
    if (a->c != NULL) VectorPrint(a->c);
    printf("A1:\n");
    if (a->A1 != NULL) MatrixPrint(a->A1);
    printf("b1: ");
    if (a->b1 != NULL) VectorPrint(a->b1);
    printf("A2:\n");
    if (a->A2 != NULL) MatrixPrint(a->A2);
    printf("b2: ");
    if (a->b2 != NULL) VectorPrint(a->b2);
}
