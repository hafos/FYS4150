#include <iomanip> // Writing to file
#include <string> // Turn to string
#include <sstream> // String formatting

#include "Lattice.hpp"
#include "timecheck.hpp"
#include "markov_chain_mc.hpp"

using namespace arma;

int main()
{
  // Test Lattice and MCMC:
  int L = 2;
  vec T = {0.1, 1.5, 3.0, 4.5, 7.0, 9.5, 30.0}; // [J/kB]

  // Sample spin configurations:
  ivec n_cycles = {100, 200, 500}; // n monte carlo cycles
  int n_steps = 100; // steps in each monte carlo cycle
  bool Ordered = 0;

  // Open outputfile
  std::ofstream ofile;
  ofile.open("simple_expectation_value_tests.dat");
  ofile << "T   n_cycles   n_steps   epsilon   m   Cv   X" << std::endl;

  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Compute exp values for different temperatures and n_cycles/steps and write to file
  for (int j=0; j<n_cycles.n_elem; j++)
  {
    for (int i=0; i<T.n_elem; i++)
    {
      double energy;
      double magnetization;
      double heat_capacity;
      double susceptibility;
      expectation_values(L, T(i), Ordered, n_cycles(j), n_steps, energy, magnetization, heat_capacity, susceptibility);
      ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << n_cycles(j)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << n_steps
      << std::setw(width) << std::setprecision(decimals) << std::scientific << energy
      << std::setw(width) << std::setprecision(decimals) << std::scientific << magnetization
      << std::setw(width) << std::setprecision(decimals) << std::scientific << heat_capacity
      << std::setw(width) << std::setprecision(decimals) << std::scientific << susceptibility
      << std::endl;
    }
  }
  // Close output file
  ofile.close();


  //std::cout << " <eps> " << energy << "  J" << endl;
  //std::cout << " <|m|> " << magnetization << endl;
  //std::cout << " Cv " << heat_capacity << "  kB" << endl;
  //std::cout << " X " << susceptibility << "  1/J" << endl;

  //timecheck();
}
