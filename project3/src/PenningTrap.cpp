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

  //Coulomb constant
  k_e = 1.38935333*pow(10, 5);
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
  //define variables
  double q_i = particles_.at(i).charge();
  vec r_i = particles_.at(i).position;
  double q_j = particles_.at(j).charge();
  vec r_j = particles_.at(j).position;

  double r = norm(r_i - r_j, 2); //absolute distance

  vec F_p = k_e * q_i * q_j * (r_i - r_j)/pow(r, 3);

  return F_p;
}

// The total force on particle_i from the external fields
vec PenningTrap::total_force_external(int i)
{
  Particle p = particles_.at(i);
  double q = p.charge();
  vec r = p.position;
  vec v = p.velocity;

  vec F_e = q*external_E_field(r);
  vec F_b = q*cross(v, external_B_field(r));

  return F_e + F_b;
}

// The total force on particle_i from the other particles
vec PenningTrap::total_force_particles(int i)
{
  mat F_p_all(3, particles_.size());
  for (int j = 0; j < particles_.size(); j++)
  {
    if (j == i)
    {
      F_p_all.col(j).zeros();
    }
    else
    {
      F_p_all.col(j) = force_particle(i, j);
    }
  }
  return sum(F_p_all, 1);
}

// The total force on particle_i from both external fields and other particles
vec PenningTrap::total_force(int i)
{
  return total_force_external(i) + total_force_particles(i);
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
  // Temporary matrices to store r,v while looping through particles
  mat r_new(3, particles_.size(), fill::zeros);
  mat v_new(3, particles_.size(), fill::zeros);

  // Calculate r_new and v_new for each particle
  for (int i = 0; i < particles_.size(); i++)
  {
    vec f = total_force(i);
    vec r = particles_.at(i).position;
    vec v = particles_.at(i).velocity;
    double m = particles_.at(i).mass();

    r_new.insert_cols(i, r + dt*v);
    v_new.insert_cols(i, v + dt*f/m);
  }

  // Then update all particles at the same time with new r and v
  for (int i = 0; i < particles_.size(); i++)
  {
    particles_.at(i).position = r_new.col(i);
    particles_.at(i).velocity = v_new.col(i);
  }
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{
  // This method currently stores a lot of extra stuff...
  // Could probably be made more effective
  // These matrices are necessary in order to update all particles simultaneously:
  mat vk1(3, particles_.size(), fill::zeros);
  mat vk2(3, particles_.size(), fill::zeros);
  mat vk3(3, particles_.size(), fill::zeros);
  mat vk4(3, particles_.size(), fill::zeros);
  mat rk1(3, particles_.size(), fill::zeros);
  mat rk2(3, particles_.size(), fill::zeros);
  mat rk3(3, particles_.size(), fill::zeros);
  mat rk4(3, particles_.size(), fill::zeros);
  // These vectors are needed to store old and intermediate states:
  std::vector<Particle> particles_old = particles_;
  std::vector<Particle> particles_new = particles_;

  // First loop: k1
  for (int i = 0; i < particles_.size(); i++)
    {
      // Get particle properties at old state
      vec f = total_force(i);
      vec r = particles_.at(i).position;
      vec v = particles_.at(i).velocity;
      double m = particles_.at(i).mass();

      //f has 3 rows and 1 column -> want insert_cols to a col with 3 rows
      vk1.insert_cols(i, dt*f/m); // slope of v at old * dt
      rk1.insert_cols(i, dt*v); // slope of r at old * dt

      // Set new particle properties at k1 (start from old)
      particles_new.at(i).position = r+0.5*rk1.col(i);
      particles_new.at(i).velocity = v+0.5*vk1.col(i);
    }
  // Update particles to k1
  particles_ = particles_new;

  // Second loop: k2
  for (int i = 0; i < particles_.size(); i++)
    {
      // Get particle properties at k1 state
      vec f = total_force(i);
      vec r = particles_.at(i).position;
      vec v = particles_.at(i).velocity;
      double m = particles_.at(i).mass();

      vk2.insert_cols(i, dt*f/m); // slope of velocity at k1 * dt
      rk2.insert_cols(i, dt*v); // slope of r at k1 * dt

      // Set new particle properties at k2 (start from old)
      r = particles_old.at(i).position;
      v = particles_old.at(i).velocity;
      particles_new.at(i).position = r+0.5*rk2.col(i);
      particles_new.at(i).velocity = v+0.5*vk2.col(i);
    }
  // Update particles to k2
  particles_ = particles_new;

  // Third loop: k3
  for (int i = 0; i < particles_.size(); i++)
    {
      // Get particle properties at k2 state
      vec f = total_force(i);
      vec r = particles_.at(i).position;
      vec v = particles_.at(i).velocity;
      double m = particles_.at(i).mass();

      vk3.insert_cols(i, dt*f/m); // slope of velocity at k2 * dt
      rk3.insert_cols(i, dt*v); // slope of r at k2 * dt

      // Set new particle properties at k3 (start from old)
      r = particles_old.at(i).position;
      v = particles_old.at(i).velocity;
      particles_new.at(i).position = r+rk3.col(i);
      particles_new.at(i).velocity = v+vk3.col(i);
    }
  // Update particles to k3
  particles_ = particles_new;

  // Fourth loop: k4
  for (int i = 0; i < particles_.size(); i++)
    {
      // Get particle properties at k3 state
      vec f = total_force(i);
      vec r = particles_.at(i).position;
      vec v = particles_.at(i).velocity;
      double m = particles_.at(i).mass();

      vk4.insert_cols(i, dt*f/m); // slope of velocity at k3 * dt
      rk4.insert_cols(i, dt*v); // slope of r at k3 * dt
    }

  // Final loop: Update particles with RK4 method
  for (int i = 0; i < particles_.size(); i++)
    {
      vec r = particles_old.at(i).position;
      vec v = particles_old.at(i).velocity;

      vec r_new = r + (rk1.col(i) + 2.*rk2.col(i) + 2.*rk3.col(i) + rk4.col(i))/6.;
      vec v_new = v + (vk1.col(i) + 2.*vk2.col(i) + 2.*vk3.col(i) + vk4.col(i))/6.;

      particles_.at(i).position = r_new;
      particles_.at(i).velocity = v_new;
    }
}
