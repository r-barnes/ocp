// Copyright (c) 2009,2010,2011,2012 Brian C. Fabien
// All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the license.txt file.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "runtime.h"
#include "vector.h"
#include "matrix.h"
#include "adiff.h"

int ADiff::number_adiff_variables = 0;
int ADiff::compute_gradient = YES;
int ADiff::compute_hessian = YES;

void ADiff::set_number_adiff_variables(int num) {
    number_adiff_variables = num;
}

void ADiff::set_compute_gradient(int yes_no) {
    compute_gradient = yes_no;
}

void ADiff::set_compute_hessian(int yes_no) {
    compute_hessian = yes_no;
}

ADiff::ADiff(const ADiff &c) {
    value = c.value;
    grad = new double[number_adiff_variables];

    hess = new double * [number_adiff_variables];
    for (register int i = 0; i < number_adiff_variables; i++)
        hess[i] = new double[number_adiff_variables];

    for (register int i = 0; i < number_adiff_variables; i++) {
        grad[i] = c.grad[i];
        for (register int j = 0; j < number_adiff_variables; j++)
            hess[i][j] = c.hess[i][j];
    }
}

ADiff::ADiff(double y, int var_num) {
    if (number_adiff_variables == 0) {
        //RuntimeWarning("number_adiff_variables = 0");
        value = 0;
        grad = NULL;
        hess = NULL;
        return;
    }

    if (var_num > number_adiff_variables)
        RuntimeError("variable number greater than number_adiff_variables");

    value = y;
    grad = new double[number_adiff_variables];

    hess = new double * [number_adiff_variables];
    for (register int i = 0; i < number_adiff_variables; i++)
        hess[i] = new double[number_adiff_variables];

    for (register int i = 0; i < number_adiff_variables; i++) {
        grad[i] = 0.0;
        for (register int j = 0; j < number_adiff_variables; j++)
            hess[i][j] = 0.0;
    }

    if (var_num > 0)
        grad[var_num - 1] = 1.0;
}

double ADiff::get_value(void) {
    return (value);
}

void ADiff::set_value(double x) {
    value = x;
}

void ADiff::set_gradient_value(int index, double x) {
    if ((index >= 0) && (index < number_adiff_variables))
        grad[index] = x;
    else
        RuntimeError("index greater than number_adiff_variables");
}

void ADiff::get_gradient(Vector g) {
    int i;
    for (i = 0; i < number_adiff_variables; i++) {
        g->e[i] = (*this).grad[i];
    }
}

void ADiff::get_hessian(Matrix H) {
    int i, j;
    for (i = 0; i < number_adiff_variables; i++) {
        for (j = 0; j < number_adiff_variables; j++) {
            H->e[i][j] = (*this).hess[i][j];
        }
    }
}

ADiff& ADiff::operator=(const ADiff &right) {
    int i, j;

    if ((this) == (&right))
        return (*this);

    (*this).value = right.value;

    for (i = 0; i < number_adiff_variables; i++) {

        (*this).grad[i] = right.grad[i];

        if (compute_hessian == YES)
            for (j = 0; j < number_adiff_variables; j++)
                (*this).hess[i][j] = right.hess[i][j];
    }
    return (*this);
}

ADiff& ADiff::operator=(double right) {
    int i, j;

    (*this).value = right;

    for (i = 0; i < number_adiff_variables; i++) {

        (*this).grad[i] = 0;

        if (compute_hessian == YES)
            for (j = 0; j < number_adiff_variables; j++)
                (*this).hess[i][j] = 0.0;
    }
    return (*this);
}

const ADiff operator-(const ADiff &left) {
    // unary minus
    int i, j;
    ADiff ans;

    ans.value = -left.value;
    for (i = 0; i < left.number_adiff_variables; i++) {

        ans.grad[i] = -left.grad[i];

        if (left.compute_hessian == YES)
            for (j = 0; j < left.number_adiff_variables; j++)
                ans.hess[i][j] = -left.hess[i][j];
    }
    return (ans);
}

const ADiff operator+(const ADiff &left) {
    // unary plus
    int i, j;
    ADiff ans;

    ans.value = left.value;
    for (i = 0; i < left.number_adiff_variables; i++) {

        ans.grad[i] = left.grad[i];

        if (left.compute_hessian == YES)
            for (j = 0; j < left.number_adiff_variables; j++)
                ans.hess[i][j] = left.hess[i][j];
    }
    return (ans);
}

const ADiff ADiff::operator+(const ADiff &right) {

    int i, j;
    ADiff ans;

    ans.value = this->value + right.value;

    for (j = 0; j < number_adiff_variables; j++) {

        ans.grad[j] = this->grad[j] + right.grad[j];

        if (compute_hessian == YES)
            for (i = 0; i < number_adiff_variables; i++)
                ans.hess[j][i] = this->hess[j][i] + right.hess[j][i];
    }
    return (ans);
}

const ADiff ADiff::operator+(double right) {
    int i, j;
    ADiff ans;

    ans.value = (*this).value + right;
    for (i = 0; i < number_adiff_variables; i++) {

        ans.grad[i] = (*this).grad[i];

        if (compute_hessian == YES)
            for (j = 0; j < number_adiff_variables; j++)
                ans.hess[i][j] = (*this).hess[i][j];
    }
    return (ans);
}

const ADiff operator+(double left, const ADiff &right) {
    int i, j;
    ADiff ans;

    ans.value = left + right.value;

    for (i = 0; i < right.number_adiff_variables; i++) {

        ans.grad[i] = right.grad[i];

        if (right.compute_hessian == YES)
            for (j = 0; j < right.number_adiff_variables; j++)
                ans.hess[i][j] = right.hess[i][j];
    }
    return (ans);
}

const ADiff operator+(const ADiff &left, const ADiff &right) {
    int i, j;
    ADiff ans;

    ans.value = left.value + right.value;

    for (i = 0; i < left.number_adiff_variables; i++) {

        ans.grad[i] = left.grad[i] + right.grad[i];

        if (left.compute_hessian == YES)
            for (j = 0; j < left.number_adiff_variables; j++)
                ans.hess[i][j] = left.hess[i][j] + right.hess[i][j];
    }
    return (ans);
}

const ADiff ADiff::operator-(const ADiff &right) {
    int i, j;
    ADiff ans;

    ans.value = (*this).value - right.value;

    for (i = 0; i < number_adiff_variables; i++) {
        ans.grad[i] = (*this).grad[i] - right.grad[i];

        if (compute_hessian == YES)
            for (j = 0; j < number_adiff_variables; j++)
                ans.hess[i][j] = (*this).hess[i][j]
                    - right.hess[i][j];
    }
    return (ans);
}

const ADiff ADiff::operator-(double right) {
    int i, j;
    ADiff ans;

    ans.value = (*this).value - right;

    for (i = 0; i < number_adiff_variables; i++) {

        ans.grad[i] = (*this).grad[i];

        if (compute_hessian == YES)
            for (j = 0; j < number_adiff_variables; j++)
                ans.hess[i][j] = (*this).hess[i][j];
    }
    return (ans);
}

const ADiff operator-(double left, const ADiff &right) {
    int i, j;
    ADiff ans;

    ans.value = left - right.value;

    for (i = 0; i < right.number_adiff_variables; i++) {

        ans.grad[i] = -right.grad[i];

        if (right.compute_hessian == YES)
            for (j = 0; j < right.number_adiff_variables; j++)
                ans.hess[i][j] = -right.hess[i][j];
    }
    return (ans);
}

const ADiff operator-(const ADiff &left, const ADiff &right) {
    int i, j;
    ADiff ans;

    ans.value = left.value - right.value;

    for (i = 0; i < left.number_adiff_variables; i++) {

        ans.grad[i] = left.grad[i] - right.grad[i];

        if (left.compute_hessian == YES)
            for (j = 0; j < left.number_adiff_variables; j++)
                ans.hess[i][j] = left.hess[i][j] - right.hess[i][j];
    }
    return (ans);
}

const ADiff ADiff::operator*(const ADiff &right) {
    int i, j;
    ADiff ans;

    ans.value = (*this).value * right.value;

    for (i = 0; i < (*this).number_adiff_variables; i++) {
        ans.grad[i] = right.value * (*this).grad[i]
                + (*this).value * right.grad[i];

        if (compute_hessian == YES)
            for (j = 0; j < (*this).number_adiff_variables; j++)
                ans.hess[i][j] = (*this).value * right.hess[i][j]
                    + right.value * (*this).hess[i][j]
                    + right.grad[i]*(*this).grad[j]
                    + (*this).grad[i] * right.grad[j];
    }
    return (ans);
}

const ADiff ADiff::operator*(double right) {
    int i, j;
    ADiff ans;

    ans.value = (*this).value*right;

    for (i = 0; i < (*this).number_adiff_variables; i++) {
        ans.grad[i] = right * (*this).grad[i];

        if (compute_hessian == YES)
            for (j = 0; j < (*this).number_adiff_variables; j++)
                ans.hess[i][j] = right * (*this).hess[i][j];
    }
    return (ans);
}

const ADiff operator*(double left, const ADiff &right) {
    int i, j;
    ADiff ans;

    ans.value = left * right.value;
    for (i = 0; i < right.number_adiff_variables; i++) {
        ans.grad[i] = left * right.grad[i];

        if (right.compute_hessian == YES)
            for (j = 0; j < right.number_adiff_variables; j++)
                ans.hess[i][j] = left * right.hess[i][j];
    }
    return (ans);
}

const ADiff operator*(const ADiff &left, const ADiff &right) {
    int i, j;
    ADiff ans;

    ans.value = left.value * right.value;

    for (i = 0; i < left.number_adiff_variables; i++) {
        ans.grad[i] = right.value * left.grad[i]
                + left.value * right.grad[i];

        if (left.compute_hessian == YES)
            for (j = 0; j < left.number_adiff_variables; j++)
                ans.hess[i][j] = left.value * right.hess[i][j]
                    + left.hess[i][j] * right.value
                    + left.grad[i] * right.grad[j]
                    + left.grad[j] * right.grad[i];
    }
    return (ans);
}

const ADiff ADiff::operator/(const ADiff &b) {
    int i, j;
    ADiff ans;
    double fac0, fac1, fac2, fac3;

    if (fabs(b.value) <= DBL_EPSILON)
        RuntimeWarning("divide by zero: ADiff operator/");

    fac0 = (1.0 / b.value);
    fac1 = ((*this).value / (b.value * b.value));
    fac2 = (1.0 / (b.value * b.value));
    fac3 = (2.0 * (*this).value / (b.value * b.value * b.value));

    ans.value = fac0 * (*this).value;

    for (i = 0; i < (*this).number_adiff_variables; i++) {

        ans.grad[i] = fac0 * (*this).grad[i] - fac1 * b.grad[i];

        if (compute_hessian == YES)
            for (j = 0; j < (*this).number_adiff_variables; j++)
                ans.hess[i][j] = fac0 * (*this).hess[i][j]
                    - fac1 * b.hess[i][j]
                    - fac2 * b.grad[j]*(*this).grad[i]
                    - fac2 * (*this).grad[j] * b.grad[i]
                    + fac3 * b.grad[i] * b.grad[j];
    }
    return (ans);
}

const ADiff ADiff::operator/(double b) {
    int i, j;
    ADiff ans;
    double fac0;

    if (fabs(b) <= DBL_EPSILON)
        RuntimeWarning("divide by zero: ADiff operator/");

    fac0 = (1.0 / b);

    ans.value = fac0 * (*this).value;

    for (i = 0; i < (*this).number_adiff_variables; i++) {

        ans.grad[i] = fac0 * (*this).grad[i];

        if (compute_hessian == YES)
            for (j = 0; j < (*this).number_adiff_variables; j++)
                ans.hess[i][j] = fac0 * (*this).hess[i][j];
    }
    return (ans);
}

const ADiff operator/(double a, const ADiff &b) {
    int i, j;
    ADiff ans;
    double fac0, fac1, fac3;

    if (fabs(b.value) <= DBL_EPSILON)
        RuntimeWarning("divide by zero: ADiff operator/");

    fac0 = (1.0 / b.value);
    fac1 = (a / (b.value * b.value));
    fac3 = (2.0 * a / (b.value * b.value * b.value));

    ans.value = fac0*a;

    for (i = 0; i < b.number_adiff_variables; i++) {

        ans.grad[i] = -fac1 * b.grad[i];

        if (b.compute_hessian == YES)
            for (j = 0; j < b.number_adiff_variables; j++)
                ans.hess[i][j] = -fac1 * b.hess[i][j]
                    + fac3 * b.grad[i] * b.grad[j];
    }
    return (ans);
}

const ADiff operator/(const ADiff &a, const ADiff &b) {
    int i, j;
    ADiff ans;
    double fac0, fac1, fac2, fac3;

    if (fabs(b.value) <= DBL_EPSILON)
        RuntimeWarning("divide by zero: ADiff operator/");

    fac0 = (1.0 / b.value);
    fac1 = (a.value / (b.value * b.value));
    fac2 = (1.0 / (b.value * b.value));
    fac3 = (2.0 * a.value / (b.value * b.value * b.value));

    ans.value = fac0 * a.value;

    for (i = 0; i < a.number_adiff_variables; i++) {

        ans.grad[i] = fac0 * a.grad[i] - fac1 * b.grad[i];

        if (a.compute_hessian == YES)
            for (j = 0; j < a.number_adiff_variables; j++)
                ans.hess[i][j] = fac0 * a.hess[i][j] - fac1 * b.hess[i][j]
                    - fac2 * b.grad[j] * a.grad[i]
                    - fac2 * a.grad[j] * b.grad[i]
                    + fac3 * b.grad[i] * b.grad[j];
    }
    return (ans);
}

const ADiff Deriv(double g, double g_u, double g_uu, const ADiff &x) {
    int i, j;
    ADiff ans;

    ans.value = g;
    for (i = 0; i < x.number_adiff_variables; i++) {
        ans.grad[i] = g_u * x.grad[i];

        if (x.compute_hessian == YES)
            for (j = 0; j < x.number_adiff_variables; j++)
                ans.hess[i][j] = g_u * x.hess[i][j] + g_uu * x.grad[i] * x.grad[j];
    }
    return (ans);
}

const ADiff sin(const ADiff &u) {

    double g, g_u, g_uu;

    g = sin(u.value);
    g_u = cos(u.value);
    g_uu = -g;

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff cos(const ADiff &u) {

    double g, g_u, g_uu;

    g = cos(u.value);
    g_u = -sin(u.value);
    g_uu = -g;

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff tan(const ADiff &u) {

    double g, g_u, g_uu, tmp;

    g = tan(u.value);
    tmp = 1.0 / cos(u.value);
    g_u = tmp*tmp; // sec^2 u
    g_uu = 2.0 * g*g_u; // 2 tan u * sec^2 u

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff sec(const ADiff &u) {

    double g, g_u, g_uu, tmp;

    tmp = 1.0 / cos(u.value);
    g = tmp; // sec u
    g_u = tmp * tan(u.value); // sec u * tan u
    g_uu = tmp * (1.0 + 2.0 * tan(u.value) * tan(u.value)); // sec u *(1 + 2*tan^2 u)

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff exp(const ADiff &u) {
    double g, g_u, g_uu;

    g = exp(u.value);
    g_u = g;
    g_uu = g;

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff log(const ADiff &u) {

    double g, g_u, g_uu;

    g = log(u.value);
    g_u = 1.0 / u.value;
    g_uu = -1.0 / (u.value * u.value);

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff sqrt(const ADiff &u) {

    double g, g_u, g_uu;

    g = sqrt(u.value);
    g_u = 0.5 / g;
    g_uu = -0.25 / (g * g * g);

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff sinh(const ADiff &u) {

    double g, g_u, g_uu;

    g = sinh(u.value);
    g_u = cosh(u.value);
    g_uu = g;

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff cosh(const ADiff &u) {

    double g, g_u, g_uu;

    g = cosh(u.value);
    g_u = sinh(u.value);
    g_uu = g;

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff abs(const ADiff &u) {
    double g, g_u, g_uu;

    if (fabs(u.value) <= DBL_EPSILON)
        RuntimeWarning("ADiff.abs(): divide by zero");

    g = fabs(u.value);
    g_u = u.value / g;
    g_uu = 0.0;

    return (Deriv(g, g_u, g_uu, u));
}

const ADiff pow(const ADiff &left, const ADiff &right) {
    ADiff tmp1, tmp2, tmp3;

    tmp1 = log(left);
    tmp2 = tmp1*right;
    tmp3 = exp(tmp2);
    return (tmp3);
}

const ADiff pow(const ADiff &left, double right) {

    double g, g_u, g_uu;

    g = pow(left.value, right);
    g_u = right * pow(left.value, right - 1.0);
    g_uu = (right - 1.0) * right * pow(left.value, right - 2.0);

    return (Deriv(g, g_u, g_uu, left));
}

const ADiff pow(double left, const ADiff &right) {
    double g, g_u, g_uu, tmp;

    g = pow(left, right.value);
    tmp = log(left);
    g_u = g*tmp;
    g_uu = g_u*tmp;

    return (Deriv(g, g_u, g_uu, right));
}

double delta(const ADiff &t, double a, double b, double c) {

    if (t.value < a)
        return (b);

    return (c);
}
