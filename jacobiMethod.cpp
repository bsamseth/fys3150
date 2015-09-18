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
	  *l = i;
	  *k = j;
	  tmp = abs(A[i][j]);
      }
    }
  }
  if(k>l){
    int tmp_l;
    tmp_l = *l;
    *l = *k;
    *k = tmp_l;
  }
}

void jacobiMethod(double** A, int n, double* lambda, double eps){
  int l,k;
  double tau,t,c,s;
  double tmp_Akk, tmp_All, tmp_Akl, tmp_Ail, tmp_Aik;
  while(off_diagonal_sum(A, n) > eps){
    find_indexes_of_max(A, n, &l, &k);
    tau = (A[l][l]-A[k][k]) / (2*A[k][l]);
    t = min(abs(-tau+sqrt(1+tau*tau)), abs(-tau-sqrt(1+tau*tau)));
    c = 1/sqrt(1+tau*tau);
    s = t*c;
    
    tmp_Akk = A[k][k];
    tmp_All = A[l][l];
    tmp_Akl = A[k][l];

    A[k][k] = tmp_Akk*c*c - 2*tmp_Akl*c*s + tmp_All*s*s;
    A[l][l] = tmp_All*c*c + 2*tmp_Akl*c*s + tmp_Akk*s*s;
    A[k][l] = A[l][k] = 0;
    for(int i = 0; i < n; i++){
      if(i!=k && i!=l){
	tmp_Aik = A[i][k];
	tmp_Ail = A[i][l];
	A[i][k] = tmp_Aik*c - tmp_Ail*s;
	A[i][l] = tmp_Ail*c + tmp_Aik*s;
      }
    }
      
  }
  for(int i = 0; i < n; i++){
    lambda[i] = A[i][i];
  }
}

int main(int argc, char ** argv){

  double** A = new double*[2];
  A[0] = new double[2];
  A[1] = new double[2];
    

  A[0][0] = 2;
  A[0][1] = -4;
  A[1][0] = -1;
  A[1][1] = -1;  
  double* lambda = new double[2];
  double eps = 1e-13;
  int l,k;
  jacobiMethod(A, 2, lambda, eps);
  cout << lambda[0] << "  " << lambda[1] << endl;
}
