// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   cholesky.h
 * Author: fabien
 *
 * Created on July 9, 2010, 6:25 PM
 */

#ifndef _CHOLESKY_H
#define	_CHOLESKY_H

#ifdef	__cplusplus
extern "C" {
#endif

int CholeskyFactor(Matrix A, Matrix L);
int CholeskySolveLinearSystem(Matrix A, Vector b, Vector x);

#ifdef	__cplusplus
}
#endif

#endif	/* _CHOLESKY_H */
