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


void tridiagonal( double b, double c, double * d, double * v ) {
    /* Function for solving the eq. Av=d in the case that A is tridiagonal,
       with a constant value on the diagonals. Specialy the subdiagonal must be equal to -1.

       b - value for the main diagonal
       c - value for the superdiagonal
       d - a vector with the function values
       v - souluiton vector. Will be filled when function returns
       N - dimension of the matrix/vectors

       Using the Thomas Algorithm, with modifications.
    */
    int i;
    double * c_prime = new double[N];

    /* The loop below will depend on c_prime[i-1] and d[i-1] to be already changed,
       so c_prime[0] and d[0] will be set manually at first */
    
    c_prime[1] = c / b;
    d[1] = d[1] / b;

    double m;
    for (i = 2; i < N-1; i++) {
	m = 1.0 / (b + c_prime[i-1]);
	c_prime[i] = c * m;
	d[i] = (d[i] + d[i-1]) * m;
    }

    v[N-2] = d[N-2];
    for (i = N - 2; i-- > 1; ) {
	v[i] = d[i] - c_prime[i] * v[i+1];
    }
    delete[] c_prime;
    delete[] d;
}



void writeToFile(string filename, double * v) {
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
    double b = 2;
    double c = -1;

    /* Make x-array. x = {0, ... , 1} with xi = i*h */
    double * x = new double[N];
    for (int i = 0; i < N; i++) {
        x[i] = i * h;
    }

    
    /* Make d-array. di = h^2*f(xi). */
    double * d = new double[N];
    for (int i = 0; i < N; i++) {
        d[i] = h*h * f(x[i]);
    }

    double * v = new double[N]; // init. solution array
    tridiagonal(b,c,d,v);
    
    stringstream outputname;
    outputname << "data/v_solve_N_" << N << ".dat";
    writeToFile(outputname.str(), v);
    cout << "Data writen to " << outputname.str() << endl;
    return 0;
}

