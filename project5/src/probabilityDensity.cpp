#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>
#include "Body.h"
#include "Universe.h"
#include "../../lib/lib.h"

using namespace arma;
using namespace std;

int main(int argc, char *argv[]) {

  bool useVerlet = true;
  
  if (argc < 5) {
    cout << "Usage: ./main.x h tmax N check_for_ejections [verlet/rk4]" << endl;
    exit(1);
  }

  double h = atof(argv[1]);
  double t_max = atof(argv[2]);
  int N = atoi(argv[3]);
  bool check_for_ejections = (bool) atoi(argv[4]);
  
  if (argc > 5) {
    string arg5 (argv[5]);
    if (arg5 == "rk4")
      useVerlet = false;
  }

  Universe mysystem;

  long idum = -1;
  
  
  const double R0 = 20;
  const double M0 = 1;
  const double t_crunch = 1;
  
  double wanted_totalmass = 1000;
  for (int i = 0; i < N; i++) {
    double mass = abs(10*M0 + M0 * gaussian_deviate(&idum));
    
    double vx, vy, vz;
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
  double factor = wanted_totalmass/totalmass;
  
  for (int i = 0; i < N; i++) {
    mysystem.all_bodies[i].mass *= factor;
  }
  
  double rho0 = wanted_totalmass / (4.0 * M_PI * pow(R0, 3) / 3);
  
  mysystem.G = 3*M_PI / (32 * t_crunch*t_crunch * rho0);
  mysystem.R0 = R0;


  t_max = t_max*t_crunch;
  
  double E0 = mysystem.energy();
  if (useVerlet){
    mysystem.solve_Verlet(h, t_max, check_for_ejections);
  }
  else
    mysystem.solve_RK4(h, t_max, check_for_ejections);
  
  double E1 = mysystem.energy();

  cout << "E0 = " << E0 << ", E1 = " << E1 << ", Delta E = " << E1-E0 << ", rel. err. = " << (E1-E0)/E0 <<  endl;
}
