#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>
#include "Body.h"
#include "Universe.h"
#include "../../lib/lib.h"

using namespace arma;
using namespace std;

int main() {
  Universe mysystem;

  long idum = -1;
  
  const int N = 100;
  const double R0 = 20;
  const double M0 = 1;
  const double t_crunch = 1;
  
  for (int i = 0; i < N; i++) {
    double mass = abs(10*M0 + M0 * gaussian_deviate(&idum));
    
    int vx, vy, vz;
    vx = vy = vz = 0;

    double phi   = 2*M_PI * ran2(&idum);
    double r     = R0 * pow(ran2(&idum), 1/3.);
    double theta = acos( 1 - 2*ran2(&idum) );

    double x = r * sin(theta)*cos(phi);
    double y = r * sin(theta)*sin(phi);
    double z = r * cos(theta);
    Body b = Body(mass, x, y, z, vx, vy, vz);
    mysystem.add_body(b);
  }

  double totalmass = 0;
  for (int i = 0; i < N; i++) {
    totalmass += mysystem.all_bodies[i].mass;
  }
  double rho0 = totalmass / (4.0 * M_PI * pow(R0, 3) / 3);
  
  mysystem.G = 3*M_PI / (32 * t_crunch*t_crunch * rho0);



  double h = 0.0011;
  double t_max = 0.9*t_crunch;

  
  double E0 = mysystem.energy();
  mysystem.solve_RK4(h, t_max);
  double E1 = mysystem.energy();

  cout << "E0 = " << E0 << ", Delta E = " << E1-E0 << ", rel. err. = " << (E1-E0)/E0 <<  endl;
    
}
