#include "Particle.hpp"

//constructor
Particle::Particle(double q_in, double m_in, vec r_in, vec v_in)
{
  charge = q_in;
  mass = m_in;
  position = r_in;
  velocity = v_in;
}

//methods that return the four variables
double Particle::charge()
{
  return charge;
}
double Particle::mass()
{
  return mass;
}
vec Particle::position()
{
  return position;
}
vec Particle::velocity()
{
  return velocity;
}
