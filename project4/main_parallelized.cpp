#include <iomanip> // Writing to file
#include <string> // Turn to string
#include <sstream> // String formatting
#include <ctime> // Time measurements
#include "omp.h"  // OpenMP header

#include "Lattice.hpp"
#include "markov_chain_mc.hpp"

using namespace arma;

int main()
{
  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Problem 8
  // Run final estimations with parallelized code
  ivec L = {40};//, 60, 80 100};
  vec T = linspace(2.1, 2.4, 2);
  bool Ordered = 0; // Less random fluctuations in computation time?
  int n_burnin = 1000;
  int n_cycles = 1000;

  for (int j=0; j<L.n_elem; j++)
  {
    int N = L(i)*L(i);

    vec energy(T.n_elem);
    vec magnetization(T.n_elem);
    vec heat_capacity(T.n_elem);
    vec susceptibility(T.n_elem);

    // Parallelized loop over temperature:

    #pragma omp parallel for
    for (int i=0; i<T.n_elem; i++)
    {
      expectation_values(L(i), T(i), Ordered, n_cycles, n_burnin, energy(i),
                          magnetization(i), heat_capacity(i), susceptibility(i));
    }
    // End of parallelization

    // Open outputfile
    std::ostringstream lstring;
    lstring << L(i);
    std::string filename = "results_L" + lstring.str() + ".dat"
    ofile.open(filename);
    ofile << "n_cycles: "
    << std::setw(width) << std::setprecision(decimals) << std::scientific << n_cycles << std::endl;
    ofile << "T  epsilon   m   Cv   X" << std::endl;
    for (int i=0; i<T.n_elem; i++)
    {
      ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << energy(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << magnetization(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << heat_capacity(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << susceptibility(i)
      << std::endl;
    }
    // Close output file
    ofile.close();

  }
}
