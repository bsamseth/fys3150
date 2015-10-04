#include <iostream>
#include <cstdlib>
#include <cmath>
#include "jacobiMethod.h"
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>

using namespace std;

// define functions used by main
// function bodys follows below
double potentialInteract(double rho, double omega_r);
double potentialNonInteract(double rho, double omega_r);
void solve_and_write(int argc, char** argv, double (*potential)(double, double), bool interacting);


int main(int argc, char** argv) {

  if (argc < 5){
    cout << "Usage interactingElectrons.cpp interacting nstep rho_max omega_r" << endl;
    cout << "interacting = 1, non-interacting = 0" << endl;
    exit(1);
  }

  bool interacting = atoi(argv[1]);

  if (interacting) 
      solve_and_write(argc, argv, potentialInteract, true);
  else
      solve_and_write(argc, argv, potentialNonInteract, false);
}


double potentialInteract(double rho, double omega_r){
    /* Potential for the interacting case */
    return pow(omega_r*rho, 2) + 1.0/rho;
}

double potentialNonInteract(double rho, double omega_r){
    /* Potential for the non-interacting case */
    return pow(omega_r*rho, 2);
}

void solve_and_write(int argc, char** argv, double (*potential)(double, double), bool interacting) {
    /* Solve the SL with potential (*potential) with parameters as
       given in argv, and write result to data file */
  double rho_min = 0;
  double rho_max = atof(argv[3]);
  int nstep = atoi(argv[2]);
  double* omega_r = new double[argc-4];
  for(int i = 4; i < argc; i++){
    omega_r[i-4] = atof(argv[i]);
  }
  
  
  double h = (rho_max - rho_min) / nstep;
  double ** A = new double*[nstep-1];
  for (int i = 0; i < nstep-1; i++) {
    A[i] = new double[nstep-1];
    for (int j = 0; j < nstep-1; j++) {
      A[i][j] = 0;
    }
  }

  for(int s = 0; s < argc-4; s++){

    double const_off_diag = -1/(h*h);
    double const_diag     = 2/(h*h);
    for (int i = 0; i < (nstep-1) -1; i++) {
	A[i][i]   = const_diag +  (*potential)(rho_min+(i+1)*h, omega_r[s]);
	A[i][i+1] = A[i+1][i] = const_off_diag;
    }
    A[nstep-2][nstep-2] = const_diag + (*potential)(rho_min+((nstep-2)+1)*h,omega_r[s]);  
    

    double* lambdas = new double[nstep-1];
    double** R = make_identity_matrix(nstep - 1);
    jacobiMethod(A, R, nstep-1, lambdas, 1e-10);

    // write solutions to a data file
    stringstream filename;
    filename << "data/";
    if (interacting == false)
	filename << "Non";
    filename << "interactingElectrons_omegar_" << omega_r[s] << "_nstep_" << nstep << "_rhomax_" << rho_max << ".dat";
    ofstream outFile;
    outFile.open(filename.str());

    vector<double> myvec (lambdas, lambdas + nstep-1);
    sort(myvec.begin(), myvec.end());
    int * columns = new int[3];
    for(int i = 0; i < 3; i++){
	outFile << myvec[i] << endl;
      for (int j = 0; j < nstep -1; j++) {
	  if (lambdas[j] == myvec[i]) {
	      columns[i] = j;
	  }
      }
    }

    for (int k = 0; k < 3; k++ ){
	for (int i = 0; i < nstep-1; i++) {
	    outFile << R[i][columns[k]] << " ";
	}
    	outFile << endl;
    } 
    outFile.close();
  }
}
