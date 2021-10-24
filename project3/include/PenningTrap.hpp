#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <vector>
#include "Particle.hpp"

using namespace arma;

class PenningTrap
{
  private:

    // Stores all the particles in the PenningTrap
    std::vector<Particle> particles_;

    // The parameters
    double B0_;
    double c_;
    double c0_;
    double d_;
    bool Interactions_;
    double time_;
    double f_;
    double wV_;

    // Constants
    double k_e_;

  public:
    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, bool Interactions);

    // Returns the particles
    std::vector<Particle> particles_in_trap();

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // Count particles in the trap
    int count_particles();

    // Set the time-dependent potential parameters:
    void set_time_dependence(double f, double wV, double start_time);

    // Update time_ and c_ (for time dependent potential):
    void update_time(double h);

    // External electric field at point r=(x,y,z)
    vec external_E_field(vec r);

    // External magnetic field at point r=(x,y,z)
    vec external_B_field(vec r);

    // Force on particle_i from particle_j
    vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    vec total_force_external(int i);

    // The total force on particle_i from the other particles
    vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    vec total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);
};

#endif
