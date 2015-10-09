#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>

#include <omp.h>

#define ZERO 1E-10

using namespace std;
using namespace std::chrono;

void gauleg(double x1, double x2, double x[], double w[], int n);
double int_function(double x1, double x2, double y1, double y2, double z1, double z2);


int main(int argc, char** argv) {
    if (argc < 4) {
	cout << "Usage: " << argv[0] << " a b N" << endl;
	exit(1);
    }

    double a = atof(argv[1]);
    double b = atof(argv[2]);
    int    N = atoi(argv[3]);


    high_resolution_clock::time_point t1 = high_resolution_clock::now();

  
    // array for integration points and weights using Legendre polynomials
    double *x = new double [N];
    double *w = new double [N];
    //   set up the mesh points and weights
    gauleg(a,b,x,w, N);
    //   evaluate the integral with the Gauss-Legendre method
    //   Note that we initialize the sum
    double int_gauss = 0.;
    //   six-double loops

    
    int i,j,k,l,m,n;
#pragma omp parallel for reduction(+:int_gauss) private(i,j,k,l,m,n)
    for (i=0;i<N;i++){
	for (j = 0;j<N;j++){
	    for (k = 0;k<N;k++){
		for (l = 0;l<N;l++){
		    for (m = 0;m<N;m++){
			for (n = 0;n<N;n++){
			    int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]
				*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
			}
		    }
		}
	    }
	}
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>( t2 - t1 ).count() / 1e9;
    
    cout << "Computed = " << int_gauss << endl;
    cout << "Exact = " << 5*M_PI*M_PI/(16.0*16) << endl;
    cout << "error = " << abs(int_gauss - 5*M_PI*M_PI/(16.0*16)) << endl;
    cout << "Time spent = " << duration << endl;
    return 0;
}


double int_function(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double alpha = 2.;
    // evaluate the different terms of the exponential
    double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
    double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
    double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    // Cheating!!
    if(deno <pow(10.,-6.)) { return 0;}
    else return exp(exp1+exp2)/deno;
} // end of function to evaluate
