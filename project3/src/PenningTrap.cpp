#include <vector>
#include <armadillo>

#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  B0 = B0_in;
  V0 = V0_in;
  d = d_in;
  std::vector<Particle> particles;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
  particles.push_back(p_in);
}

// External electric field at point r=(x,y,z)
vec PenningTrap::external_E_field(vec r)
{
  //this is the negative gradient of the potential
  double x = r(0);
  double y = r(1);
  double z = r(3);

  double c = -2.*V0/pow(d, 2);

  vec E_ext = c*{x/2., y/2., -z};
  return E_ext;
}

// External magnetic field at point r=(x,y,z)
vec PenningTrap::external_B_field(vec r)
{
  //this is just B0 in the z direction
  vec B_ext = {0, 0, B0};
  return B_ext;
}

// Force on particle_i from particle_j
vec PenningTrap::force_particle(int i, int j)
{
  //placeholder
  return 0;
}

// The total force on particle_i from the external fields
vec PenningTrap::total_force_external(int i)
{
  Particle p = particles.at(i);
  q = p.charge();
  m = p.mass();
  r = p.position();
  v = p.velocity();

  omega_0 = q*B0/m;
  omega_z2 = 2*q*V0/(m*pow(d, 2));

  F_e = q*external_E_field(r);
  F_b = q*cross(v, external_B_field(r));

  return F_e + F_b;
}

// The total force on particle_i from the other particles
vec PenningTrap::total_force_particles(int i)
{
  //this is a placeholder
  vec F_p = {0, 0, 0};

  return F_p;
}

// The total force on particle_i from both external fields and other particles
vec PenningTrap::total_force(int i)
{
  return total_force_external(i) + total_force_particles(i);
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
  for (int i = 0; i < particles.size(); i++)
  {
    f = total_force(i);
    r = particles.at(i).position();
    v = particles.at(i).velocity();
    m = particles.at(i).mass();
    q = particles.at(i).charge();

    r_new = r + dt*v;
    v_new = v + dt*f/m;

    Particle p_new(q, m, r_new, v_new);

    particles.at(i) = p_new;
  }
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{
  //nothing here yet
  return 0;
}
