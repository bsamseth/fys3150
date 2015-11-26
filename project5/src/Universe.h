#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <fstream>
#include <vector>
#include <armadillo>
#include "Body.h"

const double G = 4*M_PI*M_PI; // universal constant

class Universe {
  public:

  int n_bodies = 0;
  std::vector<Body> all_bodies;
  
  Universe();
  void add_body(Body b);
  void solve_RK4(double h, double t_max);
  void solve_Verlet(double h, double t_max);
  void initialize_system_matrix(arma::mat &ma);
  double force(double x, double y, double z, double M_other);
  void derivative(arma::mat &in_matrix_y, arma::mat &out_matrix_dy, int n);
  void sum_matrix(arma::mat &result, double coeff_one, arma::mat &first, double coeff_two, arma::mat &second, int n);
  void print_position(std::ofstream &ofile);
};


#endif
