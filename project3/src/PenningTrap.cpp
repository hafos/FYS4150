#include <vector>
#include <armadillo>

#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  B = B0_in;
  V = V0_in;
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
  double z = r(2);

  double c = 2.*V/pow(d, 2);

  vec E_ext = c*{x/2., y/2., -z};
  return E_ext;
}

// External magnetic field at point r=(x,y,z)
vec PenningTrap::external_B_field(vec r)
{
  //this is just B0 in the z direction
  vec B_ext = {0, 0, B};
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
  Particle p = particles.at(i);
  q = p.charge();
  // // m = p.mass(); // not needed here
  r = p.position();
  v = p.velocity();

  // // omega_0 = q*B/m; // not needed here
  // // omega_z2 = 2*q*V/(m*pow(d, 2)); // not needed here

  F_e = q*external_E_field(r);
  F_b = q*cross(v, external_B_field(r));

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
  // // return particles; // This function returns void
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
//unfinished
{
  // Places to store stuff in local: size 3xN_particles
  std::vector<std::vector<double>> vk1;
  std::vector<std::vector<double>> vk2;
  std::vector<std::vector<double>> vk3;
  std::vector<std::vector<double>> vk4;
  std::vector<std::vector<double>> rk1;
  std::vector<std::vector<double>> rk2;
  std::vector<std::vector<double>> rk3;
  std::vector<std::vector<double>> rk4;
  std::vector<Particle> particles_old = particles; // To store old particles
  std::vector<Particle> particles_new; // To store intermediate states while updating individual particles
  for (int i = 0; i < particles.size(); i++)
  {
    // Get particle properties at old state
    f = total_force(i);
    r = particles.at(i).position();
    v = particles.at(i).velocity();
    m = particles.at(i).mass();
    q = particles.at(i).charge();

    vk1.push_back(f/m); // slope of v at old
    rk1.push_back(v); // slope of r at old

    // Set new particle properties at k1
    Particle p_new(q, m, r+0.5*rk1.at(i), v+0.5*vk1.at(i));
    particles_new.at(i) = p_new;
  }
  // Update particles to k1
  particles = particles_new;
  for (int i = 0; i < particles.size(); i++)
  {
    // Get particle properties at k1 state
    f = total_force(i);
    r = particles.at(i).position();
    v = particles.at(i).velocity();
    m = particles.at(i).mass();
    q = particles.at(i).charge();

    vk2.push_back(f/m); // slope of velocity at k1
    rk2.push_back(v); // slope of r at k1

    // Set new particle properties at k2
    r = particles_old.at(i).position();
    v = particles_old.at(i).velocity();
    m = particles_old.at(i).mass();
    q = particles_old.at(i).charge();
    Particle p_new(q, m, r+0.5*rk2.at(i), v+0.5*vk2.at(i));
    particles_new.at(i) = p_new;
  }
  // Update particles to k2
  particles = particles_new;
  for (int i = 0; i < particles.size(); i++)
  {
    // Get particle properties at k2 state
    f = total_force(i);
    r = particles.at(i).position();
    v = particles.at(i).velocity();
    m = particles.at(i).mass();
    q = particles.at(i).charge();

    vk3.push_back(f/m); // slope of velocity at k2
    rk3.push_back(v); // slope of r at k2

    // Set new particle properties at k3
    r = particles_old.at(i).position();
    v = particles_old.at(i).velocity();
    m = particles_old.at(i).mass();
    q = particles_old.at(i).charge();
    Particle p_new(q, m, r+rk3.at(i), v+vk3.at(i));
    particles_new.at(i) = p_new;
  }
  // Update particles to k3
  particles = particles_new;
  for (int i = 0; i < particles.size(); i++)
  {
    // Get particle properties at k3 state
    f = total_force(i);
    r = particles.at(i).position();
    v = particles.at(i).velocity();
    m = particles.at(i).mass();
    q = particles.at(i).charge();

    vk4.push_back(f/m); // slope of velocity at k3
    rk4.push_back(v); // slope of r at k3
  }
  // Update particles with RK4 method
  for (int i = 0; i < particles.size(); i++)
  {
    r = particles_old.at(i).position();
    v = particles_old.at(i).velocity();
    m = particles_old.at(i).mass();
    q = particles_old.at(i).charge();

    r_new = r + (1/6)*dt*(rk1.at(i) + 2*rk2.at(i) + 2*rk3.at(i) + rk4.at(i));
    v_new = v + (1/6)*dt*(vk1.at(i) + 2*vk2.at(i) + 2*vk3.at(i) + vk4.at(i));

    Particle p_new(q, m, r_new, v_new);

    particles.at(i) = p_new;
  }
}
