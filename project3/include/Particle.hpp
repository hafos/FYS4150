#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

using namespace arma;

class Particle
{
private:

  // Stores the parameters of the particle that can't be changed from outside:
  double charge_;
  double mass_;

public:
  //constructor
  Particle(double q_in, double m_in, vec r_in, vec v_in);

  // Parameters of the particle that can be changed from outside:
  vec position;
  vec velocity;

  //methods to return charge and mass:
  double charge();
  double mass();
};

#endif
