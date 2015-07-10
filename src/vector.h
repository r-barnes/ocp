// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

#ifndef _VECTOR_H
#define	_VECTOR_H

#ifdef	__cplusplus
extern "C" {
#endif
struct vector_ {
	int 	r;   // number of rows
	double  *e;  // elements
};

#define VEC(x,i) x->e[i]
#define VEC_LENGTH(x) x->r

typedef struct vector_ * Vector;

long int VectorAllocated();
Vector *VectorArrayNew(int k, int n); // an array of k vectors dimension n
void    VectorArrayDelete(Vector *a, int k);
void    VectorAxpy(Vector y, double a, Vector x);
Vector  VectorCopy(Vector v);
double  VectorDot(Vector y, Vector x);
Vector  VectorLinspace(double a, double b, int n);
void	VectorDelete(Vector v);
int		VectorFprintf(Vector v, char *filename);
double  VectorGet(Vector v, int i);
double  VectorMax(Vector y);
double  VectorMin(Vector y);
Vector 	VectorNew(int r);
double  VectorNorm(Vector y);
double  VectorNormInf(Vector y);
double  VectorNormOne(Vector y);
void	VectorPrint(Vector v);
void    VectorSet(Vector v, int i, double x);
void    VectorSetAllTo(Vector v, double x);
void    VectorSetEqual(Vector a, Vector b);
void    VectorScale(Vector v, double s);

#ifdef	__cplusplus
}
#endif

#endif	/* _VECTOR_H */
