#ifndef FUNCTIONS_H
#define FUNCTIONS_H

const double EPS = 3.0e-14;
const double ZERO = 1E-10;
const int MAXIT = 10;

void gauleg(double x1, double x2, double x[], double w[], int n);
void gauss_laguerre(double *x, double *w, int n, double alf);
double gammln( double xx);


#endif
