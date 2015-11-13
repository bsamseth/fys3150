#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <sstream>
#include <string>
// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void initialize(int, int **, double&, double&, int, long&);
// Function to calculate E and M, but keeping the configuration
void initialize_keepconfig(int n_spins, int **spin_matrix, double& E, double& M);
// The metropolis algorithm 
void Metropolis(int, long&, int **, double&, double&, double *, int* );
// prints to file the results of the calculations  
std::string output(int, int, double, double *, int*);
//  Matrix memory allocation
//  allocate space for a matrix
void  **matrix(int, int, int);
//  free space for  a matrix
void free_matrix(void **);
// ran2 for uniform deviates, initialize with negative seed.
double ran2(long *);

#endif
