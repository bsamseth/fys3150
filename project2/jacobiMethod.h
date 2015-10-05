#ifndef JACOBI_METHOD_
#define JACOBI_METHOD_

double off_diagonal_sum(double** A, int n);
void find_indexes_of_max(double** A, int n, int* l, int* k);
double max_Aij(double** A, int n);
int jacobiMethod(double** A, double** R, int n, double* lambda, double eps);
double** make_identity_matrix(int n);

#endif
