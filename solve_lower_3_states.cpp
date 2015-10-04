#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include "jacobiMethod.h"

using namespace std;

int main(int argc, char** argv) {

    if (argc < 3) {
	cout << "Usage: a.x nstep rho_max" << endl;
	exit(0);
    }

    double tolerance = 1e-10; // for jacobiMethod, not accuracy of egenvalues
    int nstep = atoi(argv[1]);
    double rho_max = atof(argv[2]);
    double rho_min = 0.0;

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

    double * lambdas = new double[nstep-1];
    int number_of_transformations = jacobiMethod(A, make_identity_matrix(nstep-1), nstep-1, lambdas, tolerance);
    vector<double> myvec (lambdas, lambdas + nstep-1);
    sort(myvec.begin(), myvec.end());
    for (int i = 0; i < 3; i++)
    	cout << myvec[i] << endl;
    cout << "Used " << number_of_transformations << " similarity transformations." << endl;
}
