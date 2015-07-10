// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

#ifndef _QR_H
#define	_QR_H

#ifdef	__cplusplus
extern "C" {
#endif

void _householder(double **Q, double **R, double **B, int n, int m);

#ifdef	__cplusplus
}
#endif

#endif
