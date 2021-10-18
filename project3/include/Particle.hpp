#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

using namespace arma;

class Particle
{
public:
  //constructor
  Particle(double q_in, double m_in, vec r_in, vec v_in);

  //methods to return all the variables:
  double charge();
  double mass();
  vec position();
  vec velocity();
};

#endif
