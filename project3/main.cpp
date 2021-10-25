#include <iostream>
#include <iomanip> // Writing to file in nice format
#include <string> // Turn to string
#include <sstream> // String formatting
#include <armadillo> // vectors and matrices
#include <cmath>

#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;

void single_particle_test(double dt)
{
  // Set parameters
  double sim_time = 100; // [microseconds]
  double B0_in = 9.65*1e1;
  double V0_in = 9.65*1e8;
  double d_in = 1e4; // [micrometers]
  bool Inter = 0; // One particle, does not matter if on or off here

  // Initialize the system
  PenningTrap PT(B0_in, V0_in, d_in, Inter);
  PenningTrap PT2(B0_in, V0_in, d_in, Inter);
  // Add particles
  double q_in = 1;
  double m_in = 1;
  vec r_in = { 10, 0, 10 };
  vec v_in = { 0, -1, 0 };

  Particle p_in = Particle(q_in, m_in, r_in, v_in);
  PT.add_particle(p_in);
  PT2.add_particle(p_in);

  // Make test output files
  std::ostringstream dt_string;
  dt_string << std::setw(5) << std::setprecision(1) << std::scientific << dt;
  std::string filename = "test_results/test_RK4-FE_dt" + dt_string.str() + ".dat";

  std::ofstream ofile;
  ofile.open(filename);

  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Write parameters to first line of file:
  ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << B0_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << V0_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << d_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << q_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << m_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << v_in.at(1)
    << std::endl;
  ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << 0
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(0)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(1)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(2)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(0)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(1)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(2)
    << std::endl;

  // Test methods
  double t = 0;
  vec r(3);
  vec r2(3);
  bool stop_FE = 0;
  while (t <= sim_time){
    // Update:
    PT.evolve_RK4(dt);
    if (norm(PT2.particles_in_trap()[0].position) < d_in)
    { // After this point it is pointless to continue updating, it just grows
      PT2.evolve_forward_Euler(dt);
      stop_FE = 1;
    }
    t += dt;
    // Write to file:
    r = PT.particles_in_trap()[0].position;
    r2 = PT2.particles_in_trap()[0].position;
    ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << t
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r.at(0)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r.at(1)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r.at(2)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2.at(0)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2.at(1)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2.at(2)
        << std::endl;
    // Stop conditions:
    if ((norm(PT.particles_in_trap()[0].position) > d_in) and (stop_FE == 1))
    { // After this point it is pointless to continue updating, both are out of bounds
      break;
    }
    if (dt <= 0)
    {
      std::cout << "What are you doing... dt should be positive!   dt: "<< dt << std::endl;
      break;
    }
  }
  // Close file
  ofile.close();
}

void multiple_particle_test(double dt, bool Interactions)
{ // This could definitely be made prettier with some loops
  // Set parameters
  double sim_time = 100; // [microseconds]
  double B0_in = 9.65*1e1;
  double V0_in = 9.65*1e8;
  double d_in = 1e4; // [micrometers]

  // Initialize the systems
  PenningTrap PT(B0_in, V0_in, d_in, Interactions);
  // Add particles
  double q_in = 1;
  double m_in = 1;
  vec r_in = { 10, 0, 10 };
  vec v_in = { 0, -1, 0 };
  vec r_in2 = {0, 10, 0};
  vec v_in2 = { 0, 0, 0 };

  Particle p_in = Particle(q_in, m_in, r_in, v_in);
  PT.add_particle(p_in);
  Particle p_in2 = Particle(q_in, m_in, r_in2, v_in2);
  PT.add_particle(p_in2);

  // Make test output file
  std::string filename = "test_results/2particles_RK4_interactions_off.dat";
  if (Interactions==1){
    filename = "test_results/2particles_RK4_interactions_on.dat";
  }
  std::ofstream ofile;
  ofile.open(filename);
  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Test method
  double t = 0;
  vec r1(3);
  vec r2(3);
  vec v1(3);
  vec v2(3);


  // Write parameters to first line of file:
  ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << B0_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << V0_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << d_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << q_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << m_in
    << std::endl;
  ofile << "t x y z vx vy vz x2 y2 z2 vx2 vy2 vz2" << endl;
  ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << 0
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(0)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(1)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in.at(2)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << v_in.at(0)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << v_in.at(1)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << v_in.at(2)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in2.at(0)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in2.at(1)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << r_in2.at(2)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << v_in2.at(0)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << v_in2.at(1)
    << std::setw(width) << std::setprecision(decimals) << std::scientific << v_in2.at(2)
    << std::endl;

  while (t <= sim_time){
    PT.evolve_RK4(dt);
    t += dt;

    r1 = PT.particles_in_trap().at(0).position;
    r2 = PT.particles_in_trap().at(1).position;
    v1 = PT.particles_in_trap().at(0).velocity;
    v2 = PT.particles_in_trap().at(1).velocity;

    ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << t
      << std::setw(width) << std::setprecision(decimals) << std::scientific << r1.at(0)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << r1.at(1)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << r1.at(2)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << v1.at(0)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << v1.at(1)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << v1.at(2)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << r2.at(0)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << r2.at(1)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << r2.at(2)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << v2.at(0)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << v2.at(1)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << v2.at(2)
      << std::endl;
    if (dt <= 0)
      {
        std::cout << "What are you doing... dt should be positive!   dt: "<< dt << std::endl;
        break;
      }
  }
  // Close file
  ofile.close();
}

void resonance_test(double f, vec wV, double dt, bool Interactions)
{
  // Test for resonance phenomena: unfinished
  // Set parameters
  double sim_time = 500; // [microseconds]
  //double sim_time = 10; //shorter time to test if function works
  double B0_in = 9.65*1e1;
  double V0_in = 9.65*1e7*0.0025;
  double d_in = 0.05*1e4; // [micrometers]

  // Initialize the system
  PenningTrap PT(B0_in, V0_in, d_in, Interactions);
  // Add particles
  double q_in = 1;
  double m_in = 40.08;
  int n = 100;
  PT.add_n_particles(n, q_in, m_in);

  //make a vector of initial conditions of all particles
  std::vector <Particle> p_initial = PT.particles_in_trap();

  // Make test output file
  std::string filename = "test_results/resonance_interactions_off_f="+std::to_string(f)+".dat";
  if (Interactions==1){
    filename = "test_results/resonance_interactions_on_f="+std::to_string(f)+".dat";
  }
  std::ofstream ofile;
  ofile.open(filename);
  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Write parameters to first line of file:
  ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << B0_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << V0_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << d_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << q_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << m_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << f
    << std::endl;
  ofile << "wV frac" << endl;

  for (int i = 0; i < wV.size(); i++)
  {
    // Set time dependence
    double t = 0; //start time
    PT.set_time_dependence(f, wV.at(i), t);

    while (t <= sim_time)
    {
      PT.evolve_RK4(dt);
      t += dt;
    }

    //count particles in trap
    int n_particles = PT.count_particles();

    //write to file
    ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << wV.at(i)
       << std::setw(width) << std::setprecision(decimals) << std::scientific << n_particles
       << std::endl;

    std::cout << n_particles << endl;

    //reset trap to initial conditions
    PT.reset_trap(p_initial);
  }

}

int main()
{
  // Test runs with single particle, using both FE and RK4
  // They save position and time to file
  // Use plot_tests.py to view the results
  double dt = 0.01;
  //single_particle_test(dt);
  //dt = 0.005;
  //single_particle_test(dt);
  //dt = 0.001;
  //single_particle_test(dt);
  //dt = 0.0005;
  //single_particle_test(dt);
  //dt = 0.0001;
  //single_particle_test(dt);

  // Test runs with two particles, using RK4
  // They save position, velocity and time to file
  // Use plot_tests.py to view the results
  // dt = 0.0001;
  // bool Interactions = 0;
  // multiple_particle_test(dt, Interactions);
  // Interactions = 1;
  // multiple_particle_test(dt, Interactions);

  // Test resonance without interactions
  dt = 0.1;
  vec f_in = {0.1, 0.4, 0.7};
  //vec f_in = {0.1};
  vec wV_in = linspace(0.2, 2.5, 115); //MHz, step size 0.02
  //vec wV_in = {0.2};
  bool Interactions = 0;
  for (double n:f_in)
  {
    std::cout << "f = " << n << endl;
    resonance_test(n, wV_in, dt, Interactions);
  }

  //test particle count
  // double sim_time = 1; // [microseconds]
  // double B0_in = 9.65*1e1;
  // double V0_in = 9.65*1e8;
  // double d_in = 1e4;
  //
  // bool Interactions = 1;
  // PenningTrap PT(B0_in, V0_in, d_in, Interactions);
  //
  // double q_in = 1;
  // double m_in = 1;
  // vec r_in = { 1, 0, 0 };
  // vec v_in = { 0, 0.1, 0 };
  // vec r_in2 = {0, 1, 0};
  // vec v_in2 = { 0, 0, 0 };
  //
  // Particle p_in = Particle(q_in, m_in, r_in, v_in);
  // PT.add_particle(p_in);
  // Particle p_in2 = Particle(q_in, m_in, r_in2, v_in2);
  // PT.add_particle(p_in2);
  //
  // double t = 0;
  // int count = 0;
  //
  // int maxiter = 10;
  // while (t <= sim_time)
  // {
  //   PT.evolve_RK4(dt);
  //   std::cout << PT.count_particles() << endl;
  //   t += dt;
  //   count +=1;
  //   if (count >= maxiter){
  //     std::cout << "Maxiter, t:" << t << endl;
  //     break;
  //   }
  // }

  //test random particles:
  // double q_in = 1;
  // double m_in = 1;
  // double B0_in = 9.65*1e1;
  // double V0_in = 9.65*1e8;
  // double d_in = 1e4;
  // int n = 10;
  // bool Interactions = 1;
  //
  // PenningTrap PT(B0_in, V0_in, d_in, Interactions);
  // PT.add_n_particles(n, q_in, m_in);
  // std::cout << PT.count_particles() << endl;
}
