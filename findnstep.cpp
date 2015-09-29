#include <iostream>
#include <cstdlib>
#include <cmath>
#include "jacobiMethod.h"

using namespace std;


bool check_eigenvalues(double * lambdas, int number_of_lambdas, double tolerance) {
    double exact0 = 3;
    double exact1 = 7;
    double exact2 = 11;

    bool s1, s2, s3;
    s1 = s2 = s3 = false;
    for(int i = 0; i < number_of_lambdas; i++){
	s1 = s1 or ( abs( lambdas[i] - exact0) < tolerance);
	s2 = s2 or ( abs( lambdas[i] - exact1) < tolerance);
	s3 = s3 or ( abs( lambdas[i] - exact2) < tolerance);
    }

    cout << "nstep-1 = " << number_of_lambdas << ":" <<endl;
    for(int i = 0; i < number_of_lambdas; i++){
	cout << lambdas[i] << " " ;
    }
    cout << endl;
    return ( s1 and s2 and s3 );
}

int main(int argc, char** argv) {

    double tolerance = 1e-4;
    
    int nstep = 200;
    
    double rho_min = 0;
    double rho_max = 100;

    if (argc < 1)
	rho_max = atof(argv[1]);

    double * lambdas = new double[1];
    do {
	delete[] lambdas; // free memory, new lambdas will be made
	double h = (rho_max - rho_min) / nstep;
	
	double ** A = new double*[nstep-1];
	for (int i = 0; i < nstep-1; i++) {
	    A[i] = new double[nstep-1];
	    for (int j = 0; j < nstep-1; j++) {
		A[i][j] = 0;
	    }
	}
		
	double const_off_diag = -1/(h*h);
	double const_diag     = 2/(h*h);
	for (int i = 0; i < (nstep-1) -1; i++) {
	    A[i][i]   = const_diag + pow(rho_min + (i+1)*h, 2);
	    A[i][i+1] = A[i+1][i] = const_off_diag;
	}
	A[nstep-2][nstep-2] = const_diag + pow(rho_min + ((nstep-2)+1)*h, 2);


	lambdas = new double[nstep-1];
	jacobiMethod(A, nstep-1, lambdas, 1e-10);

	for (int i = 0; i < nstep-1; i++) {
	    delete[] A[i];
	}
	nstep += 10;
				 
    } while ( ( check_eigenvalues(lambdas, nstep-1 -10, tolerance) == false ) and ( nstep < 100) );

    cout << "Minimum nstep = " << nstep-1 << endl;
    return 0;
}
