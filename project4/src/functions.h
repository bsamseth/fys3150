#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void initialize(int, int **, double&, double&, int, long&);
// The metropolis algorithm 
void Metropolis(int, long&, int **, double&, double&, double *, int* );
// prints to file the results of the calculations  
const char* output(int, int, double, double *, int*);
//  Matrix memory allocation
//  allocate space for a matrix
void  **matrix(int, int, int);
//  free space for  a matrix
void free_matrix(void **);
// ran2 for uniform deviates, initialize with negative seed.
double ran2(long *);

#endif
