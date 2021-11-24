#include <iomanip> // Writing to file
#include <string> // Turn to string
#include <sstream> // String formatting
#include <ctime> // Time measurements
#include <chrono>
#include "omp.h"  // OpenMP header

#include "Lattice.hpp"
#include "markov_chain_mc.hpp"

using namespace arma;

int main()
{
  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Problem 7
  // Check for time improvements
  int L = 100;
  vec T = {2.0, 2.1, 2.2, 2.3, 2.5};
  int N = L*L;
  bool Ordered = 1; // Less random fluctuations in computation time?
  int n_burnin = 1000;
  int n_cycles = 1000;

  vec energy(T.n_elem);
  vec magnetization(T.n_elem);
  vec heat_capacity(T.n_elem);
  vec susceptibility(T.n_elem);




  time_t time_start = std::time(NULL);
  auto t1 = std::chrono::high_resolution_clock::now();

  // Loop over temperature, will be paralellized later :
  for (int i=0; i<T.n_elem; i++)
  {
    expectation_values(L, T(i), Ordered, n_cycles, n_burnin, energy(i),
                        magnetization(i), heat_capacity(i), susceptibility(i));
  }
  // End of loop

  auto t2 = std::chrono::high_resolution_clock::now();
  time_t time_end = std::time(NULL);
  std::cout << "Not Parallelized: " << endl;
  std::cout << "Time used (time) " << difftime(time_end, time_start) << " seconds " <<endl;
  double timechrono = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()*0.001;
  std::cout << "Time used (chrono) " << timechrono << " seconds" << endl;

  // Open outputfile
  std::ofstream ofile;
  ofile.open("timing_parallelization_nothing.dat");
  ofile << "timeused [s]: "
  << std::setw(width) << std::setprecision(decimals) << std::scientific << timechrono << std::endl;
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



  // Repeat but paralellized:

  time_start = std::time(NULL);
  t1 = std::chrono::high_resolution_clock::now();

  // Parallelized loop over temperature:

  #pragma omp parallel for
  for (int i=0; i<T.n_elem; i++)
  {
    expectation_values(L, T(i), Ordered, n_cycles, n_burnin, energy(i),
                        magnetization(i), heat_capacity(i), susceptibility(i));
  }
  // End of parallelization

  t2 = std::chrono::high_resolution_clock::now();
  time_end = std::time(NULL);
  std::cout << "Parallelized: " << endl;
  std::cout << "Time used (time) " << difftime(time_end, time_start) << " seconds " <<endl;
  timechrono = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()*0.001;
  std::cout << "Time used (chrono) " << timechrono << " seconds" << endl;

  // Open outputfile
  ofile.open("timing_parallelization_temperatures.dat");
  ofile << "timeused [s]: "
  << std::setw(width) << std::setprecision(decimals) << std::scientific << timechrono << std::endl;
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
