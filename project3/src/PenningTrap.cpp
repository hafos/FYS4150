#include <vector>
#include <armadillo>

#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  B_ = B0_in;
  V_ = V0_in;
  d_ = d_in;
  //std::vector<Particle> particles_; // Defined in .hpp
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
  particles_.push_back(p_in);
}

// Return the particle objects
std::vector<Particle> PenningTrap::particles_in_trap()
{
  return particles_;
}

// External electric field at point r=(x,y,z)
vec PenningTrap::external_E_field(vec r)
{
  //this is the negative gradient of the potential
  double x = r(0);
  double y = r(1);
  double z = r(2);

  double c = 2.*V_/pow(d_, 2);

  vec E_ext = {c*x/2., c*y/2., -c*z};
  return E_ext;
}

// External magnetic field at point r=(x,y,z)
vec PenningTrap::external_B_field(vec r)
{
  //this is just B0 in the z direction
  vec B_ext = {0, 0, B_};
  return B_ext;
}

// Force on particle_i from particle_j
vec PenningTrap::force_particle(int i, int j)
{
  //placeholder
  vec F_p = {0, 0, 0};
  return F_p;
}

// The total force on particle_i from the external fields
vec PenningTrap::total_force_external(int i)
{
  Particle p = particles_.at(i);
  double q = p.charge();
  // // m = p.mass(); // not needed here
  vec r = p.position();
  vec v = p.velocity();

  // // omega_0 = q*B_/m; // not needed here
  // // omega_z2 = 2*q*V_/(m*pow(d, 2)); // not needed here

  vec F_e = q*external_E_field(r);
  vec F_b = q*cross(v, external_B_field(r));

  return F_e + F_b;
}

// The total force on particle_i from the other particles
vec PenningTrap::total_force_particles(int i)
{
  //this is a placeholder
  vec F_p_all = {0, 0, 0};

  return F_p_all;
}

// The total force on particle_i from both external fields and other particles
vec PenningTrap::total_force(int i)
{
  return total_force_external(i) + total_force_particles(i);
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
  for (int i = 0; i < particles_.size(); i++)
  {
    vec f = total_force(i);
    vec r = particles_.at(i).position();
    vec v = particles_.at(i).velocity();
    double m = particles_.at(i).mass();
    double q = particles_.at(i).charge();

    vec r_new = r + dt*v;
    vec v_new = v + dt*f/m;

    Particle p_new(q, m, r_new, v_new);

    particles_.at(i) = p_new;
  }
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
//unfinished
{
  //
  mat vk1(particles_.size(),3);
  // Note: rest of my attempt is located in RK4_contents_note.txt at the moment
  // (It did not allow me to compile...)

}
