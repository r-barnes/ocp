// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   equations.h
 * Author: fabien
 *
 * Created on July 9, 2010, 6:24 PM
 */

#ifndef _EQUATIONS_H
#define	_EQUATIONS_H

#ifdef	__cplusplus
extern "C" {
#endif

struct equations_ {
    void (*F)(Vector x, Vector f); // compute f(x)
    void (*DF)(Vector x, Matrix df); // compute the Jacobian df = df/dx
    // private
    Vector x1, F0, F1;
};

typedef struct equations_ * Equations;

Equations EquationsNew();
void EquationsDelete(Equations eqns);
void EquationsDF(Equations eqns, Vector x0, Matrix Df);

#ifdef	__cplusplus
}
#endif

#endif	/* _EQUATIONS_H */
