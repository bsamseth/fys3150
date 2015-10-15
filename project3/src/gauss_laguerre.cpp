#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include <omp.h>

#include "functions.h"

using namespace std;
using namespace std::chrono;

double int_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2);

int main(int argc, char** argv) {
    if (argc < 2) {
	cout << "Usage: " << argv[0] << " N" << endl;
	exit(1);
    }

    int    N = atoi(argv[1]);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

  
    // array for integration points and weights using Legendre polynomials
    double *x_theta = new double [N];
    double *w_theta = new double [N];
    double *x_phi = new double [N];
    double *w_phi = new double [N];
    double *x_r = new double [N];
    double *w_r = new double [N];

    double alf = 2;   // for lagurre
    double alpha = 2; // from wavefunc.
    
    //   set up the mesh points and weights
    gauleg(0,M_PI, x_theta, w_theta, N/2);
    gauleg(0,2*M_PI, x_phi, w_phi, N/2);
    gauss_laguerre(x_r,w_r,N,alf);
    //   evaluate the integral with the Gauss-Legendre method
    //   Note that we initialize the sum
    double int_gauss = 0.;
    //   six-double loops
    
    int i,j,k,l,m,n;

#pragma omp parallel for reduction(+:int_gauss) private(i,j,k,l,m,n)
    for (i=0;i<N;i++){
	for (j = 0;j<N;j++){
	    for (k = 0;k<N/2;k++){
		for (l = 0;l<N/2;l++){
		    for (m = 0;m<N/2;m++){
			for (n = 0;n<N/2;n++){
			    int_gauss+=w_r[i]*w_r[j]*w_theta[k]*w_theta[l]*w_phi[m]*w_phi[n]
				*int_function(x_r[i],x_r[j],x_theta[k],x_theta[l],x_phi[m],x_phi[n]);
			}
			
		    }
		    
		}
	    }
	}
    }

    int_gauss /= pow(2*alpha,4)*2*alpha;
	
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>( t2 - t1 ).count() / 1e9;
    
    cout << "Computed = " << int_gauss << endl;
    cout << "Exact = " << 5*M_PI*M_PI/(16.0*16) << endl;
    cout << "error = " << abs(int_gauss - 5*M_PI*M_PI/(16.0*16)) << endl;
    cout << "Time spent = " << duration << endl;
    return 0;
}

double int_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2) {
    double cosBeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    double tmp = r1*r1 + r2*r2 - 2*r1*r2*cosBeta;
    if ( abs( tmp ) <= ZERO )
	return 0;
    return (sin(theta1)*sin(theta2))/ sqrt(tmp);
}
