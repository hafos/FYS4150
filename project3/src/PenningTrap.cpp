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
  mat vk1(3, particles_.size(), fill::zeros);
  mat vk2(3, particles_.size(), fill::zeros);
  mat vk3(3, particles_.size(), fill::zeros);
  mat vk4(3, particles_.size(), fill::zeros);
  mat rk1(3, particles_.size(), fill::zeros);
  mat rk2(3, particles_.size(), fill::zeros);
  mat rk3(3, particles_.size(), fill::zeros);
  mat rk4(3, particles_.size(), fill::zeros);
  std::vector<Particle> particles_old = particles_; // To store old particles
  std::vector<Particle> particles_new; // To store intermediate states while updating individual particles
  for (int i = 0; i < particles_.size(); i++)
    {
      // Get particle properties at old state
      vec f = total_force(i);
      vec r = particles_.at(i).position();
      vec v = particles_.at(i).velocity();
      double m = particles_.at(i).mass();
      double q = particles_.at(i).charge();
      //std::cout << f << endl << v << endl;
      //std::cout << vk4.n_cols << vk4.n_rows << endl;
      //f has 3 rows and 1 column -> want insert_cols to a col with 3 rows
      //

      vk1.insert_cols(i, dt*f/m); // slope of v at old * dt
      rk1.insert_cols(i, dt*v); // slope of r at old * dt

      // Set new particle properties at k1
      //std::cout << "Passed" << endl;
      Particle p_new(q, m, r+0.5*rk1.col(i), v+0.5*vk1.col(i));
      particles_new.push_back(p_new);
    }
  // Update particles to k1
  particles_ = particles_new;
  for (int i = 0; i < particles_.size(); i++)
    {
      // Get particle properties at k1 state
      vec f = total_force(i);
      vec r = particles_.at(i).position();
      vec v = particles_.at(i).velocity();
      double m = particles_.at(i).mass();
      double q = particles_.at(i).charge();

      vk2.insert_cols(i, dt*f/m); // slope of velocity at k1
      rk2.insert_cols(i, dt*v); // slope of r at k1

      // Set new particle properties at k2
      r = particles_old.at(i).position();
      v = particles_old.at(i).velocity();
      m = particles_old.at(i).mass();
      q = particles_old.at(i).charge();
      Particle p_new(q, m, r+0.5*rk2.col(i), v+0.5*vk2.col(i));
      particles_new.at(i) = p_new;
    }
  // Update particles to k2
  particles_ = particles_new;
  for (int i = 0; i < particles_.size(); i++)
    {
      // Get particle properties at k2 state
      vec f = total_force(i);
      vec r = particles_.at(i).position();
      vec v = particles_.at(i).velocity();
      double m = particles_.at(i).mass();
      double q = particles_.at(i).charge();

      vk3.insert_cols(i, dt*f/m); // slope of velocity at k2
      rk3.insert_cols(i, dt*v); // slope of r at k2

      // Set new particle properties at k3
      r = particles_old.at(i).position();
      v = particles_old.at(i).velocity();
      m = particles_old.at(i).mass();
      q = particles_old.at(i).charge();
      Particle p_new(q, m, r+rk3.col(i), v+vk3.col(i));
      particles_new.at(i) = p_new;
    }
  // Update particles to k3
  particles_ = particles_new;
  for (int i = 0; i < particles_.size(); i++)
    {
      // Get particle properties at k3 state
      vec f = total_force(i);
      vec r = particles_.at(i).position();
      vec v = particles_.at(i).velocity();
      double m = particles_.at(i).mass();
      double q = particles_.at(i).charge();

      vk4.insert_cols(i, dt*f/m); // slope of velocity at k3
      rk4.insert_cols(i, dt*v); // slope of r at k3
    }
  // Update particles with RK4 method
  for (int i = 0; i < particles_.size(); i++)
    {
      vec r = particles_old.at(i).position();
      vec v = particles_old.at(i).velocity();
      double m = particles_old.at(i).mass();
      double q = particles_old.at(i).charge();

      vec r_new = r + (rk1.col(i) + 2.*rk2.col(i) + 2.*rk3.col(i) + rk4.col(i))/6.;
      vec v_new = v + (vk1.col(i) + 2.*vk2.col(i) + 2.*vk3.col(i) + vk4.col(i))/6.;

      Particle p_new(q, m, r_new, v_new);

      particles_.at(i) = p_new;
    }
}
