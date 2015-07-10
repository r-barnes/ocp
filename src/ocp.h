// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   ocp.h
 * Author: fabien
 *
 * Created on July 11, 2010, 5:46 PM
 */

#ifndef _OCP_H
#define	_OCP_H

#ifdef	__cplusplus
extern "C" {
#endif

struct ocp_ {
    // public members
    int n_states;       // Number of states > 0  
    int n_controls;     // Number of controls >= 0
    int n_parameters;   // Number of parameters >= 0 
    int n_initial;      // Number of initial conditions >= 0 
    int n_terminal;     // Number of terminal conditions >= 0 
    int n_inequality;   // Number of inequality constraints >= 0 
    void (*initial_constraints)(Vector y, Vector p, Vector Gamma);  // Computer Gamma(y,p)
    double (*terminal_constraints)(Vector y, Vector p, Vector Psi); // Compute phi, Psi
    void (*inequality_constraints)(Vector y, Vector u, Vector p, double t, Vector d); // Compute d
    double (*differential_equations)(Vector y, Vector u, Vector p, double t, Vector F); // Compute L, F
    void (*Dinitial_constraints)(Vector y, Vector p, Matrix Gammay, Matrix Gammap); // derivatives of Gamma
    void (*Dterminal_constraints)(Vector y, Vector p, Vector phiy, Vector phip, Matrix Psiy, Matrix Psip); // derivatives of phi and Psi
    void (*Dinequality_constraints)(Vector y, Vector u, Vector p, double t, Matrix dy, Matrix du, Matrix dp); // derivatives of d
    void (*Ddifferential_equations)(Vector y, Vector u, Vector p, double t, Matrix Fy, Matrix Fu, Matrix Fp, Vector Ly, Vector Lu, Vector Lp); // derivatives of L and F

    // private
    
    double _t;
    Vector _y, _u, _p, _yp, _yup;

    Equations _G;
    Matrix _dgamma;

    Equations _phiPsi;
    Matrix _dphiPsi;

    Equations _fL;
    Matrix _dfL;

    Equations _d;
    Matrix _dd;

    int _memory_initilized;
};

typedef struct ocp_ * OCP;

OCP  OCPNew();
void OCPDelete(OCP o);
void OCPFDGamma(OCP o, Vector y, Vector p, Matrix Gammay, Matrix Gammap);
void OCPFDphiPsi(OCP o, Vector y, Vector p, Vector phiy, Vector phip, Matrix Psiy, Matrix Psip);
void OCPFDfL(OCP o, Vector y, Vector u, Vector p, double t, Matrix fy, Matrix fu, Matrix fp, Vector Ly, Vector Lu, Vector Lp);
void OCPFDd(OCP o, Vector y, Vector u, Vector p, double t, Matrix dy, Matrix du, Matrix dp);
int  OCPCheck(OCP o);

// FIX: These do not belong here, should be part of the solver
int  OCPWriteData(OCP o, FILE *file_id, Vector T, Matrix Y, Matrix U, Vector P, int control_type, double tol);
int  OCPReadData(OCP o, FILE *file_id, Vector T, Matrix Y, Matrix U, Vector P);

#ifdef	__cplusplus
}
#endif

#endif	/* _OCP_H */
