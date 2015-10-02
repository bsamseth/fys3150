#include <iostream>
#include <cstdlib>
#include <cmath>
#include "jacobiMethod.h"
#include <vector>
#include <algorithm>

using namespace std;

double potential(double rho, double omega_r){
    return pow(omega_r*rho, 2) + 1.0/rho;
    }


int main(int argc, char** argv) {

  if (argc < 4){
    cout << "Usage interactingElectrons.cpp nstep rho_max omega_r" << endl;
    exit(1);
  }

  double rho_min = 0;
  double rho_max = atof(argv[2]);
  int nstep = atoi(argv[1]);
  double* omega_r = new double[argc-3];
  for(int i = 3; i < argc; i++){
    omega_r[i-3] = atof(argv[i]);
  }
  
  
  double h = (rho_max - rho_min) / nstep;
  double ** A = new double*[nstep-1];
  for (int i = 0; i < nstep-1; i++) {
    A[i] = new double[nstep-1];
    for (int j = 0; j < nstep-1; j++) {
      A[i][j] = 0;
    }
  }

  for(int s = 0; s < argc-3; s++){

    double const_off_diag = -1/(h*h);
    double const_diag     = 2/(h*h);
    for (int i = 0; i < (nstep-1) -1; i++) {
      A[i][i]   = const_diag +  potential(rho_min+(i+1)*h, omega_r[s]);
      A[i][i+1] = A[i+1][i] = const_off_diag;
    }
    A[nstep-2][nstep-2] = const_diag + potential(rho_min+((nstep-2)+1)*h,omega_r[s]);  

  
    double* lambdas = new double[nstep-1];
    double** R = make_identity_matrix(nstep - 1);
    jacobiMethod(A, R, nstep-1, lambdas, 1e-10);
  
    vector<double> myvec (lambdas, lambdas + nstep-1);
    sort(myvec.begin(), myvec.end());
    int * columns = new int[3];
    for(int i = 0; i < 3; i++){
	cout << "eigenvalue " << i << " = " << myvec[i] << endl;
      for (int j = 0; j < nstep -1; j++) {
	  if (lambdas[j] == myvec[i]) {
	      columns[i] = j;
	  }
      }
    }

    for (int k = 0; k < 3; k++ ){
	cout << "Eigenvector " << k << " = (";
	for (int i = 0; i < nstep-1; i++) {
	    cout << R[i][columns[k]] << ", ";
	}
	cout << ")" << endl;
    }
    
  }
  
  
  


  for (int i = 0; i < nstep-1; i++) {
    delete[] A[i];
  }
			     
}
