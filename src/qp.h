// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   qp.h
 * Author: fabien
 *
 * Created on July 9, 2010, 6:27 PM
 */

#ifndef _QP_H
#define	_QP_H

#ifdef	__cplusplus
extern "C" {
#endif

struct qp_ {
    Matrix H, A1, A2;
    Vector c, b1, b2;
};

typedef struct qp_ * QP;

QP QPNew();
void QPDelete(QP q);
void QPcheck(QP q, int *a);
void QPPrint(QP a);

#ifdef	__cplusplus
}
#endif

#endif	/* _QP_H */
