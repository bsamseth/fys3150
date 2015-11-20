#include "Body.h"
#include "Universe.h"
#include <armadillo>
#include <cstdio>

using namespace arma;
using namespace std;

Universe::Universe() {}

void Universe::add_body(Body b) {
  n_bodies++;
  all_bodies.push_back(b);
}


void Universe::solve_RK4(double h, double t_max) {
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
  sprintf(filename, "Planet_position_RK4_%f_%f.dat", h, t_max);
  ofstream ofile (filename);
  
  ofile.precision(5);
    

  while(t<t_max){
    // cout << "position of body1 = ";
    // for (int i = 0; i < 3; i++) {
    //   cout << all_bodies[0].position[i] << " ";
    // }
    // cout << endl;
    
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
      for(int i=0; i<3; i++){
	all_bodies[j].position[i] = y_i(j,i);
	all_bodies[j].velocity[i] = y_i(j,i+3);
      }
      
    }
    
    print_position(ofile);

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





double Universe::force(double x, double y, double z, double M_other){
  double force=0;
  double distance=0;

  distance = x*x + y*y + z*z;

  force = G*M_other/pow(distance, 1.5);

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



void Universe::print_position(std::ofstream& ofile) {
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
