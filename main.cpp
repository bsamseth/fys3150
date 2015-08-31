#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

int N;

inline double f(double x) {
    return 100 * exp(-10*x);
}


void tridiagonal(double * a, double * b, double * c, double * d, double * v, int N) {
    /* Function for solving the eq. Av=d in the case that A is tridiagonal,
       consisting of the column vectors a, b and c.
       
       a - subdiagonal vector
       b - main diagonal vector
       c - superdiagonal vector
       d - coulmn vector equal to Av
       v - vector that will contain soulution when function returns
       N - the length of v

       Using the Thomas Algorithm, with programmatic inspiration taken from the
       Wikipedia article on said algorithm on August 31st. 2015.
    */

    int i;
    
    /* The loop below will depend on c[i-1] and d[i-1] to be already changed,
       so c[0] and d[0] will be set manually at first */
    c[0] = c[0] / b[0];
    d[0] = d[0] / b[0];

    double m;
    for (i = 1; i < N; i++) {
	m = 1.0 / (b[i] - a[i] * c[i-1]);
	c[i] = c[i] * m;
	d[i] = (d[i] - a[i] * d[i-1]) * m;
    }

    v[N-1] = d[N-1];
    for (i = N - 1; i-- > 0; ) {
	v[i] = d[i] - c[i] * v[i+1];
    }
}



void writeToFile(string filename,int N, double * v) {
    ofstream outFile;
    outFile.open(filename);
    outFile << N << endl;
    for (int i = 0; i < N; i++) {
        outFile << v[i] << endl;
    }
    outFile.close();
}

int main(int argc, char ** argv) {

    if (argc < 1) {
	cout << "Usage: " << argv[0] << " dimensionN" << endl;
	exit(1);
    }

    N = atoi(argv[1]);
    
    double h = 1.0 / (N-1);

    double * a = new double[N];
    double * b = new double[N];
    double * c = new double[N];
    double * d = new double[N];

    /* first, populate the constants a, b, c with the proper values */
    for (int i = 0; i < N; i++) {
        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
    }

    /* Make x-array. x = {0, ... , 1} with xi = i*h */
    double * x = new double[N];
    for (int i = 0; i < N; i++) {
        x[i] = i * h;
    }

    /* Make d-array. di = h^2*f(xi). */
    for (int i = 0; i < N; i++) {
        d[i] = h*h * f(x[i]);
    }

    double * v = new double[N]; // init. solution array
    tridiagonal(a,b,c,d,v,N);
    stringstream outputname;
    outputname << "data/v_solve_N_" << N << ".dat";
    writeToFile(outputname.str(), N, v);
    cout << "Data writen to " << outputname.str() << endl;
    // delete [] a,b,c,d,v;
    return 0;
}

