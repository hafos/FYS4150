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

void multiple_particle_test()
{
  // Set parameters
  double sim_time = 1; // [microseconds]
  double dt = 0.001; // [microseconds]
  double B0_in = 9.65*1e1;
  double V0_in = 9.65*1e8;
  double d_in = 1e4;
  bool Interactions = 1;

  // Initialize the system
  PenningTrap PT(B0_in, V0_in, d_in, Interactions);
  // Add particles
  double q_in = 1;
  double m_in = 1;
  vec r_in = { 1, 0, 0 };
  vec v_in = { 0, 0.1, 0 };
  vec r_in2 = {0, 1, 0};
  vec v_in2 = { 0, 0, 0 };

  Particle p_in = Particle(q_in, m_in, r_in, v_in);
  PT.add_particle(p_in);
  Particle p_in2 = Particle(q_in, m_in, r_in2, v_in2);
  PT.add_particle(p_in2);

  // Make test output file
  //std::string filename = "euler_test.dat";
  //std::string filename = "RK4_test.dat";
  //std::string filename = "euler_test_2particles.dat";
  std::string filename = "RK4_test_2particles.dat";
  std::ofstream ofile;
  ofile.open(filename);
  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Test method
  double t = 0;
  int count = 0;

  //int maxiter = 1001; // Safeguard
  int maxiter = 100001; // Safeguard
  vec r1(3);
  vec r2(3);
  // Write parameters to first line of file:
  ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << B0_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << V0_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << d_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << q_in
    << std::setw(width) << std::setprecision(decimals) << std::scientific << m_in
    << endl;
  while (t <= sim_time){
    PT.evolve_RK4(dt);
    //PT.evolve_forward_Euler(dt);
    t += dt;
    count +=1;
    if (count >= maxiter){
      std::cout << "Maxiter, t:" << t << endl;
      break;
    }

    r1 = PT.particles_in_trap().at(0).position;
    r2 = PT.particles_in_trap().at(1).position;

    ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << t
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r1(0)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r1(1)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r1(2)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2(0)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2(1)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2(2)
        << std::endl;
  }
  // Close file
  ofile.close();
}

int main()
{
  // Test runs with single particle, using both FE and RK4
  // They save position and time to file
  // Use plot_tests.py to view the results
  double dt = 0.01;
  single_particle_test(dt);
  dt = 0.005;
  single_particle_test(dt);
  dt = 0.001;
  single_particle_test(dt);
  dt = 0.0005;
  single_particle_test(dt);
  dt = 0.0001;
  single_particle_test(dt);

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

  //multiple_particle_test();
}
