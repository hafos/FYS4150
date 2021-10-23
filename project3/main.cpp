#include <iostream>
#include <iomanip> // Writing to file in nice format
#include <armadillo> // vectors and matrices
#include <cmath>

#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;

int main()
{
  // Set parameters
  double sim_time = 1; // [microseconds]
  double dt = 0.0001; // [microseconds]
  double B0_in = 9.65*1e1;
  double V0_in = 9.65*1e8;
  double d_in = 1e4;

  // Initialize the system
  PenningTrap PT(B0_in, V0_in, d_in);
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
  std::string filename = "RK4_test.dat";
  //std::string filename = "euler_test_2particles.dat";
  std::ofstream ofile;
  ofile.open(filename);
  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Test method
  double t = 0;
  int count = 0;

  int maxiter_rk4 = 1001; // Safeguard
  int maxiter_euler = 100001; // Safeguard
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
    if (count >= maxiter_rk4){
    //if (count >= maxiter_euler)
      std::cout << "Maxiter, t:" << t << endl;
      break;
    }

    r1 = PT.particles_in_trap().at(0).position;
    r2 = PT.particles_in_trap().at(1).position;
    //std::cout << r << endl;

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
