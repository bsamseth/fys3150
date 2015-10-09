#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>

#include <omp.h>

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

    //   Note that we initialize the sum
    double int_MC = 0.;
    double sum_sigma = 0.;

    double invers_period = 1./RAND_MAX; // initialise the random number generator
    
    int i;
#pragma omp parallel for reduction(+:int_MC, sum_sigma) private(i)
    for (i=0;i<N;i++){
      double x1 = (double(rand())*invers_period)*(b-a)+a; 
      double x2 = (double(rand())*invers_period)*(b-a)+a; 
      double y1 = (double(rand())*invers_period)*(b-a)+a; 
      double y2 = (double(rand())*invers_period)*(b-a)+a; 
      double z1 = (double(rand())*invers_period)*(b-a)+a; 
      double z2 = (double(rand())*invers_period)*(b-a)+a; 
      double funval = int_function(x1,x2,y1,y2,z1,z2);
      int_MC += funval;
      sum_sigma += funval*funval;
    }
    double w = 1./(N);

    double standard_avvik =pow((b-a),6)*sqrt((w*sum_sigma-w*w*int_MC*int_MC));
    int_MC *= w*pow((b-a),6); //Jacobi and weights
    

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<nanoseconds>( t2 - t1 ).count() / 1e9;
    
    cout << "Computed = " << int_MC << endl;
    cout << "Standard avvik = " << standard_avvik << endl;
    cout << "Exact = " << 5*M_PI*M_PI/(16.0*16) << endl;
    cout << "error = " << abs(int_MC - 5*M_PI*M_PI/(16.0*16)) << endl;
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
