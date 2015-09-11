#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <armadillo>
using namespace std;
using namespace std::chrono;
using namespace arma;

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
    

    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    tridiagonal(b,c,d,v);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration_tridiag = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();

    mat A (N-2,N-2);
    A.fill(0);
    vec mid_diag = ones<vec>(N-2)*2;
    vec sub_diag = ones<vec>(N-3)*(-1);
    A.diag() = mid_diag;
    A.diag(1) = sub_diag;
    A.diag(-1) = sub_diag;
    mat L, U, P, y;
    lu(L, U , P, A);
    vec v1;
    vec d_copy(N-2);
    for(int i = 0; i<N-2; i++){
      d_copy[i] = f(x[i+1])*h*h;
    }
    t1 = high_resolution_clock::now();
    solve(y, P*L, d_copy);
    solve(v1, U, y);
    t2 = high_resolution_clock::now();
    
    auto duration_armaLU = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
    
    double elapsed_time_tridiag = double (duration_tridiag)/1e9;
    double elapsed_time_armaLU = double (duration_armaLU)/1e9;

    cout << N << "   " << elapsed_time_tridiag << "   " << elapsed_time_armaLU <<  endl;
    
    return 0;
}

