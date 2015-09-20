#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;




double off_diagonal_sum(double** A, int n){
  double sum = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i!=j){
      sum += A[i][j]*A[i][j];
      }
    }
  }
  assert(sum >= 0);
  return sqrt(sum);
}


void find_indexes_of_max(double** A, int n, int* l, int* k){
  *l = *k = 0;
  double tmp = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i!=j && abs(A[i][j]) > tmp){
	*k = i;  //k should be the row and l the column
	*l = j;
	tmp = abs(A[i][j]);
      }
    }
  }
}

double max_Aij(double** A, int n){
  // simple function for finding the max
  // off diagonal-element of A
  double max_ = 0.0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i != j && A[i][j]*A[i][j] > max_){
	max_ = A[i][j]*A[i][j];
      }
    }
  }
  return max_;
}

void jacobiMethod(double** A, int n, double* lambda, double eps){
  int l,k;
  double tau,t,c,s;
  double tmp_Akk, tmp_All, tmp_Akl, tmp_Ail, tmp_Aik;
  while(off_diagonal_sum(A, n) > eps){
    find_indexes_of_max(A, n, &l, &k);
    tau = (A[l][l]-A[k][k]) / (2*A[k][l]); //hva hvis A[k][l] = 0?
    if(tau > 0) {  //slik de gjorde det i lecture-notes
      t = 1.0/(tau + sqrt(1.0 + tau*tau)); 
    } else { //ekvivalent med det vi hadde?
      t = -1.0/(-tau + sqrt(1.0 + tau*tau));
    }
    c = 1.0/sqrt(1+t*t); //her hadde vi tau*tau
    s = t*c;
    
    tmp_Akk = A[k][k];
    tmp_All = A[l][l];
    tmp_Akl = A[k][l];

    A[k][k] = tmp_Akk*c*c - 2.0*tmp_Akl*c*s + tmp_All*s*s;
    A[l][l] = tmp_All*c*c + 2.0*tmp_Akl*c*s + tmp_Akk*s*s;
    A[k][l] = A[l][k] = 0;//(tmp_Akk-tmp_All)*c*s+tmp_Akl*(c*c-s*s);
    for(int i = 0; i < n; i++){
      if(i!=k && i!=l){
	tmp_Aik = A[i][k];
	tmp_Ail = A[i][l];
	A[i][k] = A[k][i] = tmp_Aik*c - tmp_Ail*s;
	A[i][l] = A[l][i] = tmp_Ail*c + tmp_Aik*s;
      }
    }
    counter++;  
  }
  for(int i = 0; i < n; i++){
    lambda[i] = A[i][i];
  }
}

int main(int argc, char ** argv){

  //some code for generating a matrix
  int n = 4;
  double** A = new double*[n];
  for(int i=0; i<n; i++){
    A[i] = new double[n];
  }

  //filling matrix with 'arbitrary' values
  int counter = 1;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      A[i][j] = counter;
      counter += 2;
    }
  }
  
  
  double* lambda = new double[n];
  double eps = 1e-15;
  int l,k;
  jacobiMethod(A, n, lambda, eps);
  for(int i = 0; i < n; i++){
    cout << lambda[i] << "  " << endl;
  }
}
