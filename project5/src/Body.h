#ifndef BODY_H
#define BODY_H

class Body {

  public:
  double position[3];
  double velocity[3];
  double mass;

  Body(double mass, double x, double y, double z, double vx, double vy, double vz);
  Body();
};



#endif
