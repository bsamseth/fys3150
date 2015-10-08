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
#define ZERO 1E-6
#define EPS 3.0e-14
#define MAXIT 10



using namespace std;
using namespace std::chrono;


double int_function(double r1, double r2, double theta1, double theta2, double phi1, double phi2);

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

double gammln(double);

//  Note that you need to call it with a given value of alpha,
// called alf here. This comes from x^{alpha} exp(-x)

void gauss_laguerre(double *x, double *w, int n, double alf)
{
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}
// end function gaulag

double gammln( double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln
#undef EPS
#undef MAXIT


void gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = M_PI;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */
 
	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /* 
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()


