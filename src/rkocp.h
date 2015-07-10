// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   rkocp.h
 * Author: fabien
 *
 * Created on July 12, 2010, 11:49 AM
 */

#ifndef _RKOCP_H
#define	_RKOCP_H

#ifdef	__cplusplus
extern "C" {
#endif

// Constant control type
#define CONSTANT 0

// Linear control type (default) 
#define LINEAR 2

// Cubic control type 
#define CUBIC 8

// Runge-Kutta methods

// Bogacki-Shampine 3(2) Runge-Kutta method 
#define BogackiShampine3_2 0

// Classical 4-th order Runge-Kutta method 
#define Classical4_3 1

// Dormand-Prince 5(4) method 
#define DormandPrince5_4 2

// Verner Robust 6(5) method
#define Verner6_5 3

// Verner Robust 7(6) method
#define Verner7_6 4

// Verner Robust 8(7) method
#define Verner8_7 5

// Verner Robust 9(8) method
#define Verner9_8 6

// Dormand-Prince 8(7) method (default)
#define DormandPrince8_7 7

// Runge-Kutta-Verner, DVERK 6(5) method 
#define DVERK6_5 8

// Runge-Kutta-Fehlberg 7(8) method
#define RungeKuttaFehlberg7_8 9

// Runge-Kutta-Merson 4("5") method
#define Merson4_5 10

// RadauIIA 3rd order implicit method
#define RadauIIA3 21

// RadauIIA 5th order implicit method
#define RadauIIA5 22

typedef struct rkocp_ * RKOCP;

struct rkocp_ {
    // Input
    
    int control_type;           // Control type: LINEAR (default)
    int runge_kutta_method;     // Runge-Kutta method: DormandPrince8_7 (default)
    enum NLPmethod nlp_method;  // NLP solver L1SQPmethod (default)
    double tolerance;           // Convergence toerance 1.0e-6 (default)
    int maximum_iterations;     // Maximum number of iterations: 5000 (default)
    int display_data;           // Show progress: NO (default)
    
    // Input/Output
    
    Vector T; // -> t0:  Time nodes: double[n_nodes] 
    Vector P; // -> p0:  Parameters: double[n_parameters]
    Matrix Y; // -> y0:  States:  double[n_nodes][n_states]
    Matrix U; // -> u0:  Controls: double[n_nodes][n_controls] 

    NLPStatistics nlp_stats;
    
    // Private members
    
    Matrix Yw; // States and cost functional: dimension[n_nodes][n_states+1]
    const char *controlTypeName;
    const char *rkMethodName;
    
    // OCP
    OCP ocp;
    int nlp_n, nlp_m, nlp_p; // NLP problem dimensions
    RKOCP prob;
    int n_nodes, n_states, n_parameters, n_controls, n_initial;
    int n_terminal, n_inequality, c_type;
    Vector yc, uc, pc;
    Vector Gamma, Psi, F, G;
    Matrix Fy, Fu, Fp, Gy, Gu, Gp;
    Matrix Gammay, Gammap, Psiy, Psip;
    Vector Ly, Lu, Lp, phiy, phip;
    
    // Integration data
    // -: f, h, g calculation
    Vector *Kstage; // stage derivatives
    Matrix *Ydata;  // stage values for the state
    Matrix *Udata;  // stage values for the controls
    Vector lte, local_error;
    double max_lte;
    int only_compute_lte;
    
    // -: sensitivity calculation, df, dh, dg
    Matrix Yx; // state sensitivity
    Matrix *Yx_stage; // sensitivity stage values
    Matrix *Kx_stage; // sensitivity stage derivatives
    
    // Runge-Kutta
    int rk_method, rk_stages, rk_order;
    Matrix rk_a, rk_w;
    Vector rk_b, rk_c, rk_e;
    
    // SQP
    Vector x0, lambda0, mu0;
    double tol;
    int max_iter;
    
    // Control interpolation coefficients
    Matrix *Ucubic;

    // Implicit Runge-Kutta methods
    Matrix *Lnode, *Unode;
    int **pLUnode;
    int integration_fatal_error;
    int column_index; // used in the integration of the sensitivity equations

    int _memory_allocated;
};

// public methods
RKOCP RKOCPNew(OCP o);
int   RKOCPSolve(RKOCP r, Vector t0, Matrix y0, Matrix u0, Vector p0);
void  RKOCPDelete(RKOCP r);

// private methods
double _RKOCPfhgE(Vector nlp_x, Vector nlp_h, Vector nlp_g);

double  _RKOCPfhg(Vector x, Vector h, Vector g);
void    _RKOCPDfhg(Vector nlp_x, Vector nlp_df, Matrix nlp_dh, Matrix nlp_dg);
int     _RKOCPcheck(RKOCP r, Vector t0, Matrix y0, Matrix u0, Vector p0);
int     _RKOCPsolveNLP(RKOCP r);
void    _RKOCPformUcubic(RKOCP r);
void    _RKOCPHermite(double *L, double t);
void    _RKOCPa2d_free_dbl(double **array, int dim1);
double  **_RKOCPa2d_allo_dbl(int dim1, int dim2);
void    _RKOCPa3d_free_dbl(double ***array, int dim1, int dim2);
double  ***_RKOCPa3d_allo_dbl(int dim1, int dim2, int dim3);
void    _getRKcoefficients0(RKOCP r);
void    _getRKcoefficients1(RKOCP r);
void    _getRKcoefficients2(RKOCP r);
void    _getRKcoefficients3(RKOCP r);
void    _getRKcoefficients4(RKOCP r);
void    _getRKcoefficients5(RKOCP r);
void    _getRKcoefficients6(RKOCP r);
void    _getRKcoefficients7(RKOCP r);
void    _getRKcoefficients8(RKOCP r);
void    _getRKcoefficients9(RKOCP r);
void    _getRKcoefficients10(RKOCP r);

void    _formInitialEstimate(RKOCP r);
void    _getSolution(RKOCP r);
int     _remesh(RKOCP r);
void    _setUdata(RKOCP r, Vector nlp_x);
void    _setUdata0(RKOCP r, Vector x);
void    _setUdata2(RKOCP r, Vector x);
void    _setUdata8(RKOCP r, Vector x);
double  _compute_LTE(RKOCP r, double *local_error, double *y0, double *y1);
void    _axpy(double *a, double x, double *y, int n);
void    _setDhGamma(RKOCP r, Vector ys, Vector ps, Matrix nlp_dh);
void    _setDhPsi(RKOCP r, Matrix Yxi, Matrix nlp_dh);
void    _setYx_stage(RKOCP r, Matrix Yx_stage_l, Matrix Yxi, int i_node, int l, double delta, Matrix *Kx_stage);
void    _Mapsx(RKOCP r, Matrix Yxi, double s, Matrix B, int i_node);
void    _setYx(RKOCP r, Matrix Yxi, Matrix *Kstagei, int i_node, double delta);
void    _setDg(RKOCP r, Matrix Yxi, Vector ys, Vector us, Vector ps, double t, int i_node, Matrix nlp_dg);
void    _addDphi(RKOCP r, Matrix Yxi, Vector nlp_df);
void    _setKx_stage(RKOCP r, Matrix Kx_stage_l,
            Matrix Yx_stage_l, int i_node, int l, Vector ys,
            Vector us, Vector ps, double t);
            
// Implicit R-K methods
void _getRKcoeffRadauIIA3(RKOCP r);
void _getRKcoeffRadauIIA5(RKOCP r);
void _formIRKstage_values(RKOCP r, Vector y0, double h, Vector *K, Matrix Y);
int _formIRKmatrixZ(RKOCP r, Vector y, Vector u, Vector p, double t, double h, Matrix L, Matrix U, int *P);
double _RKOCPfhg_implicitZ(Vector nlp_x, Vector nlp_h, Vector nlp_g);
void _RKOCPDfhg_implicitZ(Vector nlp_x, Vector nlp_df, Matrix nlp_dh, Matrix nlp_dg);
void _estimateZ(RKOCP r, Vector y0, Vector *Z0, Vector y1, double tau, Vector *Z1);
void _setZx_stage(RKOCP r, int iter, Matrix *Zx_stage);
void _setYx_stage_implicit(RKOCP r, Matrix Yx, Matrix *Zx_stage, Matrix *Yx_stage);
void _setYx_dot_stage_implicit(RKOCP r, Matrix *Zx_stage, double h, int l, Matrix Yx_dot);
void _updateZx_stage(RKOCP r, Matrix *Zx_stage, Matrix dZx);
void _setResidual_implicit(RKOCP r, int i_node, int l, Matrix Yx_dot_l, Matrix Yx_stage_l, Vector ys, Vector us, Vector ps, double t, Matrix Kx_stage_l, Matrix res);

#ifdef	__cplusplus
}
#endif

#endif	/* _RKOCP_H */

