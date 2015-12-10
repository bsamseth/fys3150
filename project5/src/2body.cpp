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

  double mass = 1;
  
  double h = 0.0001;
  double t_max = 15;
  int N = 2;
  bool check_for_ejections = (bool) 0;

  Universe mysystem;
  mysystem.G = 1;

  double vx, vy, vz;
  vx = vy = vz = 0;
  vy = -0.1;

  double x = 0;
  double y = 0;
  double z = 0;
  Body b = Body(mass, x, y, z, vx, vy, vz);
  mysystem.add_body(b);
 
  vx = vy = vz = 0;
  vy = 0.1;

  x = 5;
  y = 0;
  z = 0;
  b = Body(mass, x, y, z, vx, vy, vz);
  mysystem.add_body(b);
 

  double E0 = mysystem.energy();

  mysystem.solve_Verlet(h, t_max, true);
  
  double E1 = mysystem.energy();

  cout << "E0 = " << E0 << ", E1 = " << E1 << ", Delta E = " << E1-E0 << ", rel. err. = " << (E1-E0)/E0 <<  endl;
    
}
