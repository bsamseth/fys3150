#include "Body.h"

Body::Body(double mass, double x, double y, double z, double vx, double vy, double vz) {
  this->mass = mass; // (*this).mass = mass
  position[0] = x;
  position[1] = y;
  position[2] = z;

  velocity[0] = vx;
  velocity[1] = vy;
  velocity[2] = vz;
}

Body::Body() {}
