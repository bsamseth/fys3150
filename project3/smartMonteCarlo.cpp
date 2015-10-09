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


double int_function(double x1, double x2, double y1, double y2, double z1, double z2);


int main(int argc, char** argv) {
    if (argc < 4) {
	cout << "Usage: " << argv[0] << " a b N" << endl;
	exit(1);
    }

    double a = atof(argv[1]);
    double b = atof(argv[2]);
    int    N = atoi(argv[3]);

    double alpha = 2;
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

  
    // array for integration points and weights using Legendre polynomials

    //   Note that we initialize the sum
    double int_gauss = 0.;
    //   six-double loops

    double invers_period = 1./RAND_MAX; // initialise the random number generator
    
    int i;
#pragma omp parallel for reduction(+:int_gauss) private(i)
    for (i=0;i<N;i++){
	double u1 = -log(1 - (double(rand())*invers_period));
	double u2 = -log(1 - (double(rand())*invers_period)); 
	double theta1 = (double(rand())*invers_period)*M_PI;
	double theta2 = (double(rand())*invers_period)*M_PI;
	double phi1 = (double(rand())*invers_period)*2*M_PI; 
	double phi2 = (double(rand())*invers_period)*2*M_PI; 
	int_gauss += int_function(u1,u2,theta1,theta2,phi1,phi2);
    }
    int_gauss /= (N*pow(2.0*alpha,5));
    int_gauss *= 4*pow(M_PI,4); // jacobi 

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
    return r1*r1*r2*r2*(sin(theta1)*sin(theta2))/ sqrt(tmp);
}



