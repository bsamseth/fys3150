#include <iostream>
#include <armadillo>

#include "Body.h"
#include "Universe.h"


using namespace arma;
using namespace std;


int main()
{
Universe mysystem;

// Body body1(1, 0,0,0, 0,-0.1,0 );
// Body body2(1, 5,0,0, 0, 0.1,0 );
Body Sun(1,0,0,0,0,0,0);
Body Mercury(1.2e-7, 0.39, 0, 0,0,9.96,0);
Body Venus(2.4e-6, 0.72, 0, 0,0,7.36,0);
Body Earth(1.5e-6,1,0,0, 0, 6.26, 0);
Body Mars(3.3e-7, 1.52, 0, 0,0,5.06,0);
Body Jupiter(9.5e-4, 5.20, 0,0,0,2.75,0);
Body Saturn(2.75e-4, 9.54, 0, 0,0,2.04,0);
Body Uranus(4.4e-5, 19.19, 0, 0,0,1.43,0);
Body Neptune(5.1e-5, 30.06, 0, 0,0,1.14,0);
Body Pluto(5.6e-9, 39.53, 0, 0,0,0.99,0);


// mysystem.add_body(body1);
// mysystem.add_body(body2);

mysystem.add_body(Sun);
mysystem.add_body(Mercury);
mysystem.add_body(Venus);
mysystem.add_body(Earth);
mysystem.add_body(Mars);
mysystem.add_body(Jupiter);
mysystem.add_body(Saturn);
mysystem.add_body(Uranus);
mysystem.add_body(Neptune);
mysystem.add_body(Pluto);

int elements = mysystem.n_bodies;
cout << "number of elements = " << elements<< endl;

mysystem.solve_RK4(0.001, 248);

}
