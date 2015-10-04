#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

#include "jacobiMethod.h"

using namespace std;

// declaring functions that will be used by the main funciton
void run_all_tests();
void test_JacobiMethod_twoByTwoMatrix();
void test_JacobiMethod_threeByThreeMatrix();
void test_JacobiMethod_eigenvector();
bool check_column_vector_equal(double ** mat1, double ** mat2, int col1, int col2, int n, double eps, int phase);


// declaring variable to be updated by every test
// this is the total number of tests that run when run_all_tests() is called
int number_of_tests = 0;

int main(int argc, char **argv)
{
    run_all_tests();
    return 0;
}

void run_all_tests() {
    cout << "Running all tests." << endl;
    cout << "Inside test_JacobiMethod_twoByTwoMatrix() ...";
    test_JacobiMethod_twoByTwoMatrix();
    cout << " success" << endl;
    cout << "Inside test_JacobiMethod_threeByThreeMatrix() ...";
    test_JacobiMethod_threeByThreeMatrix();
    cout << " success" << endl;
    cout << "Inside test_JacobiMethod_eigenvector() ...";
    test_JacobiMethod_eigenvector();
    cout << " success" << endl;
    cout << "Ran " << number_of_tests << " tests. " << endl;
}

void test_JacobiMethod_twoByTwoMatrix() {
    /*
      Test that jacobiMethod returns correct eigenvalues for
      a 2x2 matrix. Test also that only 1 transformation is used.
     */
    number_of_tests++;
    // init a test matrix
    int n = 2;
    double** A = new double*[n];
    for(int i=0; i<n; i++){
	A[i] = new double[n];
    }
    /* setting a "random" symetric matrix
       A =
     
       3  7
       7  5

       lambda = 4 \pm 5 * sqrt(2) 
    */
    A[0][0] = 3; A[0][1] = 7;
    A[1][0] = 7; A[1][1] = 5;
  
    double * lambda_exact = new double[n];
    lambda_exact[0] = 4 + 5 * sqrt(2); 
    lambda_exact[1] = 4 - 5 * sqrt(2);
      
    double* lambda = new double[n];  // array for computed lambdas
    double eps = 1e-14;             
    int number_of_transformations = jacobiMethod(A, make_identity_matrix(n), n, lambda, eps); // fill lambda with solution


    /* test that values are correct
       we don't know the order of the computed lambdas
       so some extra boolean expressions are needed
    */
    bool s1, s2;
    s1 = s2 = false;
    for(int i = 0; i < n; i++){
	s1 = s1 or ( abs( lambda[i] - lambda_exact[0]) < eps);
	s2 = s2 or ( abs( lambda[i] - lambda_exact[1]) < eps);
    }
    // asserting one by one for added clarity
    assert(s1);
    assert(s2);
    // for the 2x2 case, number_of_transformations should be 1
    assert(number_of_transformations == 1);
}


void test_JacobiMethod_threeByThreeMatrix() {
    /*
      Test that jacobiMethod returns correct eigenvalues for 
      a 3x3 case.
     */
    number_of_tests++;
    // init a test matrix
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

    double * lambda_exact = new double[n];
    lambda_exact[0] = 0.0; 
    lambda_exact[1] = 4 + sqrt(46);
    lambda_exact[2] = 4 - sqrt(46);
      
    double* lambda = new double[n];  // array for computed lambdas
    double eps = 1e-14;             
    jacobiMethod(A, make_identity_matrix(n), n, lambda, eps); // fill lambda with solution


    /* test that values are correct
       we don't know the order of the computed lambdas
       so some extra boolean expressions are needed
    */
    bool s1, s2, s3;
    s1 = s2 = s3 = false;
    for(int i = 0; i < n; i++){
	s1 = s1 or ( abs( lambda[i] - lambda_exact[0]) < eps);
	s2 = s2 or ( abs( lambda[i] - lambda_exact[1]) < eps);
	s3 = s3 or ( abs( lambda[i] - lambda_exact[3]) < eps);
    }
    // asserting one by one for added clarity
    assert(s1);
    assert(s2);
    assert(s3);
}


void test_JacobiMethod_eigenvector() {
    /*
      Test that jacobiMethod returns correct eigenvectors
      for a 3x3 case.
     */
    number_of_tests++;
    // init a test matrix
    const int n = 3;
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

    /* setting up the exact eigenvectors.
       set up in a matrix M = [v1 v2 v3]
       where
       v1 = (1/15 * (1+sqrt(46)), 2/15 * (1+sqrt(46)), 1)
       v2 = (1/15 * (1-sqrt(46)), -2/15 * (-1+sqrt(46)), 1)
       v3 = (-2, 1, 0)
       These vectors need to be normalized before they can be tested.
    */
  
    double ** exact = new double*[n];
    for (int i = 0 ; i < n; i++) {
	exact[i] = new double[n];
    }
    double v1_norm = 1.0 / sqrt( pow((1.0/15 * (1+sqrt(46))), 2) + pow(2.0/15 * (1+sqrt(46)), 2) + 1);
    double v2_norm = 1.0 / sqrt( pow((1.0/15 * (1-sqrt(46))), 2) + pow(-2.0/15 * (-1+sqrt(46)), 2) + 1);
    double v3_norm = 1.0 / sqrt( 5 );
    exact[0][0] = 1.0/15 * (1+sqrt(46)) * v1_norm;
    exact[1][0] = 2.0/15 * (1+sqrt(46)) * v1_norm;
    exact[2][0] = 1 * v1_norm;
    exact[0][1] = 1.0/15 * (1-sqrt(46)) * v2_norm;
    exact[1][1] = -2.0/15 * (-1+sqrt(46)) * v2_norm;
    exact[2][1] = 1 * v2_norm;
    exact[0][2] = -2 * v3_norm;
    exact[1][2] = 1 * v3_norm;
    exact[2][2] = 0;


    double ** calculated = make_identity_matrix(n); // matrix that will store eigenvectors
    double * lambdas = new double[n];
    jacobiMethod(A, calculated, n, new double[n], 1e-14); // fill calculated with eigenvector

    
    /* We don't know the order of the eigenvectors, or the phase,
       so some extra tests has to be made.
    */
    double eps = 1e-14; // presision we require
    bool match_calculated_n[n] = { false };
    for (int calculated_n = 0; calculated_n < n; calculated_n++) {
	for (int col =0 ; col < n; col++) {
	    for (int phase = -1; phase < 2; phase += 2) {
		match_calculated_n[calculated_n] += check_column_vector_equal(exact, calculated, col, calculated_n, n, eps, phase);
	    }
	}
    }
  
    for (int i = 0 ; i < n; i++) {
	assert ( match_calculated_n[i] );
    }
}


bool check_column_vector_equal(double ** mat1, double ** mat2, int col1, int col2, int n, double eps, int phase) {
    bool success = true;
    for (int i = 0 ; i < n; i++) {
	success = success and abs( mat1[i][col1] - mat2[i][col2]*phase ) < eps;
    }
    return success;
}
