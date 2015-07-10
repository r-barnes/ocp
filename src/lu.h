// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   lu.h
 * Author: fabien
 *
 * Created on August 11, 2010, 12:30 PM
 */

#ifndef _LU_H
#define	_LU_H

#ifdef	__cplusplus
extern "C" {
#endif

int LUFactor(Matrix M, Matrix L, Matrix U, int *p);
int LUSolve(Matrix L, Matrix U, int *p, Vector y, Vector x);
int LUSolveM(Matrix L, Matrix U, int *p, Matrix Y, Matrix X);
int LUSolveLinearSystem(Matrix A, Vector b, Vector x);

#ifdef	__cplusplus
}
#endif

#endif	/* _LU_H */
