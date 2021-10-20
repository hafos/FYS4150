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
  double dt = 0.001; // [microseconds]
  double B0_in = 9.65*1e1;
  double V0_in = 9.65*1e8;
  double d_in = 1e4;

  // Initialize the system
  PenningTrap PT(B0_in, V0_in, d_in);
  // Add two particles
  double q_in = 1;
  double m_in = 1;
  vec r_in = { 0, 1, 0 };
  vec v_in = { 0.1, 0, 0 };
  vec v_in2 = { 0.01, 0, 0 };

  Particle p_in = Particle(q_in, m_in, r_in, v_in);
  PT.add_particle(p_in);
  Particle p_in2 = Particle(q_in, m_in, r_in, v_in2);
  PT.add_particle(p_in2);

  // Make test output file
  std::string filename = "euler_test.dat";
  std::ofstream ofile;
  ofile.open(filename);
  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Test euler method
  double t = 0;
  int count = 0;
  int maxiter = 1001; // Safeguard
  vec r(3);
  vec r2(3);
  while (t <= sim_time){
    PT.evolve_forward_Euler(dt);
    t += dt;
    count +=1;
    if (count >= maxiter){
      std::cout << "Maxiter, t:" << t << endl;
      break;
    }

    r = PT.particles_in_trap().at(0).position();
    r2 = PT.particles_in_trap().at(1).position();
    //std::cout << r << endl;

    ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << t
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r(0)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r(1)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r(2)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2(0)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2(1)
        << std::setw(width) << std::setprecision(decimals) << std::scientific << r2(2)
        << std::endl;
  }
  // Close file
  ofile.close();
}
