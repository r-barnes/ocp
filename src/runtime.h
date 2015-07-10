// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   runtime.h
 * Author: fabien
 *
 * Created on July 9, 2010, 6:21 PM
 */

#ifndef _RUNTIME_H
#define	_RUNTIME_H

#ifdef	__cplusplus
extern "C" {
#endif
#define BOUNDS_CHECK 0
#define YES 1
#define NO 0
#define SQRT_EPS 1.49011611938477e-08

#define NUM_CORES 8

//double fmax(double a, double b);
//double fmin(double a, double b);
double sign(double x);

double drand();
int irand(int M);
void displayTime();
void pause();
void RuntimeError(const char *string);
void RuntimeWarning(const char *string);
#ifdef	__cplusplus
}
#endif

#endif	/* _RUNTIME_H */
