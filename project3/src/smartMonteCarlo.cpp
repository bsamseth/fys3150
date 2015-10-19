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
    if (argc < 2) {
	cout << "Usage: " << argv[0] << " N" << endl;
	exit(1);
    }
    int    N = atoi(argv[1]);

    double alpha = 2;
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    //   Note that we initialize the sum
    double int_MC = 0.;
    double sum_sigma = 0;

    double invers_period = 1./RAND_MAX; // initialise the random number generator
    
    int i;
#pragma omp parallel for reduction(+:int_MC, sum_sigma) private(i)
    for (i=0;i<N;i++){
	double u1 = -log(1 - (double(rand())*invers_period));
	double u2 = -log(1 - (double(rand())*invers_period)); 
	double theta1 = (double(rand())*invers_period)*M_PI;
	double theta2 = (double(rand())*invers_period)*M_PI;
	double phi1 = (double(rand())*invers_period)*2*M_PI; 
	double phi2 = (double(rand())*invers_period)*2*M_PI; 
	double funval = int_function(u1,u2,theta1,theta2,phi1,phi2);
	int_MC += funval;
	sum_sigma += funval*funval;
    }
    double w = 1./(N);
    double standard_avvik = (4*pow(M_PI,4)/pow(2.0*alpha,5))
	*sqrt(sum_sigma*w - int_MC*int_MC*w*w);
    int_MC /= pow(2.0*alpha,5);
    int_MC *= 4*pow(M_PI,4)*w; // jacobi 

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>( t2 - t1 ).count() / 1e9;
    
    cout << "Computed = " << int_MC << endl;
    cout << "standard avvik = " << standard_avvik << endl;
    cout << "Exact = " << 5*M_PI*M_PI/(16.0*16) << endl;
    cout << "error = " << abs(int_MC - 5*M_PI*M_PI/(16.0*16)) << endl;
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



