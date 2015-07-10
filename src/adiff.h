// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

/* 
 * File:   adiff.h
 * Author: fabien
 *
 * Created on July 17, 2010, 8:31 PM
 */

#ifndef _ADIFF_H
#define	_ADIFF_H

class ADiff {
private:
    static int number_adiff_variables;
    static int compute_gradient;
    static int compute_hessian;
    friend const ADiff Deriv(double g, double g_u, double g_uu, const ADiff &x);
public:
    double value;
    double *grad;
    double **hess;

    ADiff(double y = 0, int var_num = -1);

    ~ADiff(void) {
        delete [] grad;
        for (register int i = number_adiff_variables - 1; i >= 0; i--)
            delete [] hess[i];
        delete [] hess;
    }
    ADiff(const ADiff &c);

    static void set_number_adiff_variables(int num);
    static void set_compute_gradient(int yes_no);
    static void set_compute_hessian(int yes_no);

    void get_gradient(Vector gradient);
    void get_hessian(Matrix hessian);
    double get_value(void);
    void set_value(double x);
    void set_gradient_value(int index, double x);

    ADiff & operator=(const ADiff &right);
    ADiff & operator=(double right);

    const ADiff operator+(const ADiff &right);
    const ADiff operator+(double right);
    friend const ADiff operator+(double left, const ADiff &right);
    friend const ADiff operator+(const ADiff &left, const ADiff &right);

    const ADiff operator-(const ADiff &right);
    const ADiff operator-(double right);
    friend const ADiff operator-(double left, const ADiff &right);
    friend const ADiff operator-(const ADiff &left, const ADiff &right);

    const ADiff operator*(const ADiff &right);
    const ADiff operator*(double right);
    friend const ADiff operator*(double left, const ADiff &right);
    friend const ADiff operator*(const ADiff &left, const ADiff &right);

    const ADiff operator/(const ADiff &right);
    const ADiff operator/(double right);
    friend const ADiff operator/(double left, const ADiff &right);
    friend const ADiff operator/(const ADiff &left, const ADiff &right);

    friend const ADiff operator-(const ADiff &left); // unary minus
    friend const ADiff operator+(const ADiff &left); // unary plus

    friend const ADiff abs(const ADiff &left);
    friend const ADiff sin(const ADiff &left);
    friend const ADiff cos(const ADiff &left);
    friend const ADiff tan(const ADiff &left);
    friend const ADiff sec(const ADiff &left);
    friend const ADiff exp(const ADiff &u);
    friend const ADiff log(const ADiff &u);
    friend const ADiff sinh(const ADiff &left);
    friend const ADiff cosh(const ADiff &left);
    friend const ADiff sqrt(const ADiff &u);
    friend const ADiff pow(const ADiff &left, const ADiff &right);
    friend const ADiff pow(const ADiff &left, double right);
    friend const ADiff pow(double left, const ADiff &right);

    friend double delta(const ADiff &t, double a, double b, double c);
};
#endif	/* _ADIFF_H */
