#include "Particle.hpp"

//constructor
Particle::Particle(double q_in, double m_in, vec r_in, vec v_in)
{
  charge_ = q_in;
  mass_ = m_in;
  position = r_in;
  velocity = v_in;
}

//methods that let you access the private variables
double Particle::charge()
{
  return charge_;
}
double Particle::mass()
{
  return mass_;
}
