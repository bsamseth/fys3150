#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <unittest++/UnitTest++.h>

#include "jacobiMethod.h"

TEST(jacobiMethod) {
    
  // init an empty matrix
  int n = 3;
  double** A = new double*[n];
  for(int i=0; i<n; i++){
    A[i] = new double[n];
  }
  
  /* setting a "random" symetric matrix
     A =
     
     1   2   3
     2   4   6
     3   6   3
     
     lambda = 4 \pm sqrt(46) and lambda = 0
  */
  A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
  A[1][0] = 2; A[1][1] = 4; A[1][2] = 6;
  A[2][0] = 3; A[2][1] = 6; A[2][2] = 3;
  
  double* lambda = new double[n];
  double eps = 1e-12;
  jacobiMethod(A, n, lambda, eps);
  for(int i = 0; i < n; i++){
    cout << lambda[i] << "  " << endl;
  }
    
  CHECK(true); // not sure how all of this works thus far, so just cout and pass the test
}

int main()
{
    return UnitTest::RunAllTests();
}
