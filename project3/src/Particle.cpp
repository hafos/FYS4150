#include "Particle.hpp"

//constructor
Particle::Particle(double q_in, double m_in, vec r_in, vec v_in)
{
  charge_ = q_in;
  mass_ = m_in;
  position_ = r_in;
  velocity_ = v_in;
}

//methods that return the four variables
double Particle::charge()
{
  return charge_;
}
double Particle::mass()
{
  return mass_;
}
vec Particle::position()
{
  return position_;
}
vec Particle::velocity()
{
  return velocity_;
}
