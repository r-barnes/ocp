// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "runtime.h"
#include "vector.h"

static long int _vector_allocated_ = 0;

long int VectorAllocated() {
    return (_vector_allocated_);
}

// create a r-dimensional vector
Vector VectorNew(int r) {
    Vector v;
    if (BOUNDS_CHECK)
        if (r < 1) {
            RuntimeError("VectorNew: dimension < 1");
        }
    v = (Vector) malloc(sizeof (struct vector_));
    v->r = r;
    v->e = (double *) calloc(r, sizeof (double));
    _vector_allocated_++;
    return v;
}

// create k, n-dimensional vectors
Vector *VectorArrayNew(int k, int n) {
    int i;
    Vector *a;
    if (BOUNDS_CHECK)
        if ((k < 1) || (n < 1)) {
            RuntimeError("VectorArray: invalid dimensions");
        }
    a = (Vector *) malloc(k * sizeof (Vector));
    for (i = 0; i < k; i++) {
        a[i] = VectorNew(n);
    }
    return a;
}

// copy v
Vector VectorCopy(Vector v) {
    if (v == NULL) {
        return NULL;
    }
    Vector w = VectorNew(v->r);
    int i;
    for (i = 0; i < v->r; i++) {
        w->e[i] = v->e[i];
    }
    return w;
}

// delete the array of vectors
void VectorArrayDelete(Vector *a, int k) {
    int i;
    if (a == NULL) {
        return;
    }
    for (i = 0; i < k; i++) {
        VectorDelete(a[i]);
    }
    free(a);
}

// delete the Vector data structure
void VectorDelete(Vector v) {
    if (v == NULL) {
        return;
    }
    free(v->e);
    free(v);
    _vector_allocated_--;
}

// v->e[i] = x
void VectorSet(Vector v, int i, double x) {
    if (BOUNDS_CHECK)
        if ((i < 0) || (i > v->r - 1)) {
            RuntimeError("VectorSet: index out of bounds");
        }
    v->e[i] = x;
}

// x = v->e[i]
double VectorGet(Vector v, int i) {
    if (BOUNDS_CHECK)
        if ((i < 0) || (i > v->r - 1)) {
            RuntimeError("VectorGet: index out of bounds");
        }
    return v->e[i];
}

// a = b
void VectorSetEqual(Vector a, Vector b) {
    int i;
    if (BOUNDS_CHECK)
        if (a->r != b->r) {
            RuntimeError("VectorSetEqual: invalid dimensions");
        }
    for (i = 0; i < a->r; i++) {
        a->e[i] = b->e[i];
    }
}

// v->e[i] = x, for all i
void VectorSetAllTo(Vector v, double x) {
    int i;
    for (i = 0; i < v->r; i++) {
        v->e[i] = x;
    }
}

// v->e[i] = s*v->e[i], for all i
void VectorScale(Vector v, double s) {
    int i;
    for (i = 0; i < v->r; i++) {
        v->e[i] *= s;
    }
}

// print the vector to standard output
void VectorPrint(Vector v) {
    int i;
    if (v == NULL) {
        printf("NULL\n");
        return;
    }
    for (i = 0; i < v->r; i++) {
        printf("%g ", v->e[i]);
    }
    printf("\n");
}

// print the the vector to filename
int	VectorFprintf(Vector v, char *filename) {
	FILE *file_id = NULL;
	file_id = fopen(filename, "w");
	if (file_id == NULL) {
		return -1;
	}
	if (v == NULL) {
		fprintf(file_id, "NULL\n");
		return 0;
	}
	int i;
	for (i = 0; i < v->r; i++) {
        fprintf(file_id, "%g ", v->e[i]);
    }
    fprintf(file_id, "\n");
    fclose(file_id);
    return 0;
}

// Return a vector dimension n with elements
// linearly spaced from a to b
Vector VectorLinspace(double a, double b, int n) {
    Vector t;
    double dt;
    int i;
    if (n < 2) {
        RuntimeError("VectorLinspace: invalid number of points, n < 2");
    }
    t = VectorNew(n);
    dt = (b - a) / (((double) n) - 1.0);
    t->e[0] = a;
    for (i = 1; i < n - 1; i++) {
        t->e[i] = a + ((double) i) * dt;
    }
    t->e[n - 1] = b;
    return t;
}

// y = a * x + y
void VectorAxpy(Vector y, double a, Vector x) {
    int i;
    if (BOUNDS_CHECK)
        if (y->r != x->r) {
            RuntimeError("VectorAxpy: vectors incompatible");
        }
    for (i = 0; i < y->r; i++) {
        y->e[i] += a * x->e[i];
    }
}

// s = y.dot(x)
double VectorDot(Vector y, Vector x) {
    int i;
    double s = 0.0;
    if (BOUNDS_CHECK)
        if (y->r != x->r) {
            RuntimeError("VectorDot: vectors incompatible");
        }
    for (i = 0; i < y->r; i++) {
        s += y->e[i] * x->e[i];
    }
    return s;
}

// min(y->e[i])
double VectorMin(Vector y) {
    int i;
    double x;
    x = y->e[0];
    for (i = 1; i < y->r; i++) {
        if (y->e[i] < x) {
            x = y->e[i];
        }
    }
    return x;
}

// max(y->e[i])
double VectorMax(Vector y) {
    int i;
    double x;
    x = y->e[0];
    for (i = 1; i < y->r; i++) {
        if (y->e[i] > x) {
            x = y->e[i];
        }
    }
    return x;
}

// sqrt(y.dot(y))
double VectorNorm(Vector y) {
    return sqrt(VectorDot(y, y));
}

// max(|y->e[i]|)
double VectorNormInf(Vector y) {
    int i;
    double s, x;
    x = fabs(y->e[0]);
    for (i = 1; i < y->r; i++) {
        s = fabs(y->e[i]);
        if (s > x) {
            x = s;
        }
    }
    return x;
}

// sum(|y->e[i]|)
double VectorNormOne(Vector y) {
    int i;
    double x = 0.0;
    for (i = 0; i < y->r; i++) {
        x += fabs(y->e[i]);
    }
    return x;
}
