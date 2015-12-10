#include "Body.h"
#include "Universe.h"
#include <armadillo>
#include <cstdio>
#include <cassert>

using namespace arma;
using namespace std;

Universe::Universe() {}

void Universe::add_body(Body b) {
  n_bodies++;
  all_bodies.push_back(b);
}

void Universe::solve_Verlet(double h, double t_max, bool check_for_ejections) {
  mat y_i(n_bodies,7);
  mat r_i_dt(n_bodies,7);
  mat a_dt(n_bodies,7);
  mat v_dt(n_bodies,7);
  mat v_dt_2(n_bodies,7);

  initialize_system_matrix(y_i);
  
  double t = 0;

  //print file___________________
  char *filename = new char[1000];
  sprintf(filename, "data/Nbody_position_Verlet_%f_%f_%d_%d.dat", log(h), t_max, n_bodies, check_for_ejections);
  char *filename2 = new char[1000];
  sprintf(filename2, "data/Nbody_energy_Verlet_%f_%f_%d_%d.dat", log(h), t_max, n_bodies, check_for_ejections);
  
  ofstream output (filename);
  output.precision(6);
  ofstream output2 (filename2);
  output2.precision(6);
  
  int number_of_points = t_max/h;
  int divisor1 = number_of_points/1000;
  int divisor2= number_of_points/1000;
  int counter = 0;
  while(t < t_max){

    derivative(y_i, a_dt, n_bodies);

    for(int j=0; j<n_bodies; j++){

      for(int i=0; i<3; i++){
	y_i(j,i) = y_i(j,i) + h*y_i(j,i+3) + 0.5*h*h*a_dt(j,i+3);
	v_dt_2(j,i+3) = y_i(j,i+3) + 0.5*h*a_dt(j,i+3);
      }
    }

    derivative(y_i, a_dt, n_bodies);

    for(int j=0; j<n_bodies; j++){

      for(int i=3; i<6; i++){

	y_i(j,i) = v_dt_2(j,i) + 0.5*h*a_dt(j,i);

      }
      if(!check_for_ejections or !check_ejected(j)){  //Hvis check_for_ejections er false, skjer ikke check_ejected()
	for(int i=0; i<3; i++){
	  all_bodies[j].position[i] = y_i(j,i);
	  all_bodies[j].velocity[i] = y_i(j,i+3);
	}
      }
    }
    
    print_position(output, counter%divisor1 == 0);
    print_energy(output2, counter%divisor2 == 0);

    t+=h;
    counter++;
  }
  output.close();
}

void Universe::solve_RK4(double h, double t_max, bool check_for_ejections) {
  mat y_i(n_bodies,7);
  mat y_i_temp(n_bodies,7);
  mat k1(n_bodies,7);
  mat k2(n_bodies,7);
  mat k3(n_bodies,7);
  mat k4(n_bodies,7);

  initialize_system_matrix(y_i);
  double t = 0;


  //for print file_________
  char *filename = new char[1000];
  sprintf(filename, "data/Planet_position_RK4_%f_%f_%d.dat", h, t_max,check_for_ejections);
  ofstream ofile (filename);
  
  ofile.precision(5);
    

  while(t<t_max){    
    derivative(y_i, k1, n_bodies);

    sum_matrix(y_i_temp, 1, y_i, 0.5*h, k1, n_bodies);
    derivative(y_i_temp, k2, n_bodies);
    
    sum_matrix( y_i_temp, 1,  y_i, 0.5*h,  k2, n_bodies);
    derivative( y_i_temp,  k3, n_bodies);
	
    sum_matrix( y_i_temp, 1,  y_i, h,  k3, n_bodies);
    derivative( y_i_temp,  k4, n_bodies);

    for(int j=0; j<n_bodies; j++){

      for(int i=0; i<6; i++){
	y_i(j,i) = y_i(j,i) + h*( k1(j,i) + 2*k2(j,i) + 2*k3(j,i) + k4(j,i) )/ 6.0;
      }
      
      //Syncronize position and velocity with the classes
      // Body body_tmp = all_bodies[j];
      if(!check_for_ejections or !check_ejected(j)){ //Hvis check_for_ejections er false, skjer ikke check_ejected()
	for(int i=0; i<3; i++){
	  all_bodies[j].position[i] = y_i(j,i);
	  all_bodies[j].velocity[i] = y_i(j,i+3);
	}
      }
    }
    print_position(ofile, true);

    t+=h;
    
  }

  ofile.close();
}

  
  

void Universe::initialize_system_matrix(mat &ma){
    for(unsigned i=0; i < all_bodies.size(); i++){
        Body &body_tmp = all_bodies[i];
        ma(i,6)=body_tmp.mass;

        for(int k=0; k<3;k++){
            ma(i,k)=body_tmp.position[k];
            ma(i,k+3)=body_tmp.velocity[k];
        }
    }
}

bool Universe::check_ejected(int j){
  double radius = pow(all_bodies[j].position[0],2) +
  pow(all_bodies[j].position[1],2) +
  pow(all_bodies[j].position[2],2);
  double my_energy = energy_of(j);
  if(my_energy >= 0 && radius > R0*R0){
    all_bodies[j].mass = 0;
    all_bodies[j].velocity[0] = 0;
    all_bodies[j].velocity[1] = 0;
    all_bodies[j].velocity[2] = 0;
    all_bodies[j].position[0] = 0;
    all_bodies[j].position[1] = 0;
    all_bodies[j].position[2] = -R0 -R0*(1e-3*(j+1));
    
    return true;
  }
  return false;
}

double Universe::force(double x, double y, double z, double M_other){
  double force=0;
  double distance=0;
  distance = x*x + y*y + z*z;

  force = G*M_other/pow(distance + eps, 1.5);

  return force;
}
void Universe::derivative(mat &in_matrix_y, mat &out_matrix_dy, int n){

  double accelleration_x=0,accelleration_y=0,accelleration_z=0, mod_force;
  for(int i=0; i<n; i++){

    accelleration_x=0,accelleration_y=0,accelleration_z=0;
    for(int j=0; j<n; j++){
      if(i!=j){

	mod_force = force(in_matrix_y(j,0)-in_matrix_y(i,0),
			  in_matrix_y(j,1)-in_matrix_y(i,1),
			  in_matrix_y(j,2)-in_matrix_y(i,2),
			  in_matrix_y(j,6));


	accelleration_x += mod_force*(in_matrix_y(j,0)-in_matrix_y(i,0));
	accelleration_y += mod_force*(in_matrix_y(j,1)-in_matrix_y(i,1));
	accelleration_z += mod_force*(in_matrix_y(j,2)-in_matrix_y(i,2));
      }
    }
    out_matrix_dy(i,3) = accelleration_x;
    out_matrix_dy(i,4) = accelleration_y;
    out_matrix_dy(i,5) = accelleration_z;

  }


  for(int i=0; i<n; i++){
    out_matrix_dy(i,0) = in_matrix_y(i,3); //velx
    out_matrix_dy(i,1) = in_matrix_y(i,4); //vely
    out_matrix_dy(i,2) = in_matrix_y(i,5); //velz
  }
}




void Universe::sum_matrix(mat &result, double coeff_one, mat &first,double coeff_two, mat &second, int n){
  for(int j=0; j<n; j++){
    for(int i=0; i<6; i++){
      result(j,i) = coeff_one*first(j,i) + coeff_two*second(j,i);
    }
    result(j,6) = first(j,6);
  }
}



void Universe::print_position(std::ofstream& ofile, bool do_print) {
  if(!do_print) {return;}
  int n=3;
  for(unsigned i=0; i < all_bodies.size(); i++){
    Body &body_tmp = all_bodies[i];
    // std::cout << std::scientific;
    for(int j=0; j<n;j++){
      // std::cout << body_tmp.position[j] << "   ";
      ofile << std::scientific << body_tmp.position[j] << "   ";
    }
    // std::cout << "         ";
    ofile  << "         ";
  }
  // std::cout << std::endl;
  ofile << endl;
}


void Universe::print_energy(std::ofstream& ofile, bool do_print) {
 if(do_print){
   double Ek = 0;
   for(int j = 0; j < n_bodies; j++){
     Ek += energy_kinetic_of(j);
   }
   ofile << std::scientific << energy()<< " ";
   ofile << std::scientific << Ek << endl;
 }
}


double Universe::energy_kinetic_of(int j) {
  Body me = all_bodies[j];
  return 0.5*me.mass* (pow(me.velocity[0], 2) + pow(me.velocity[1], 2) + pow(me.velocity[2], 2));
}

double Universe::energy_potential_of(int j) {
  Body you, me = all_bodies[j];
  double E = 0;
  for (int i = 0; i < n_bodies; i++) {
    if (i == j) continue;
    you = all_bodies[i];
    double distance = pow(me.position[0]-you.position[0],2) +
      pow(me.position[1]-you.position[1],2) +
      pow(me.position[2]-you.position[2],2) + eps;
    
    assert (distance > 0);
    E += - G * you.mass * me.mass / sqrt(distance);
  }
  return E;
}

double Universe::energy_of(int j) {
  return energy_kinetic_of(j) + energy_potential_of(j);
}


double Universe::energy() {
  Body me, you;
  double E = 0;
  int N = n_bodies;
  for (int i = 0; i < N; i++) {
    me = all_bodies[i];
    E += energy_kinetic_of(i);
    for (int j = i+1; j < N; j++) {
      you = all_bodies[j];
      double distance = pow(me.position[0]-you.position[0],2) +
	pow(me.position[1]-you.position[1],2) +
	pow(me.position[2]-you.position[2],2) + eps;
      
      assert (distance > 0);
      E += - G * you.mass * me.mass / sqrt(distance);
    }
  }
  return E;
}
