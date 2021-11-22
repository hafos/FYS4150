#include <iomanip> // Writing to file
#include <string> // Turn to string
#include <sstream> // String formatting

#include "Lattice.hpp"
#include "markov_chain_mc.hpp"

using namespace arma;

int main()
{
  // Problem 4:
  int L = 2;
  vec T = {0.8, 1.0, 1.5, 2.5, 4.5, 7.0, 10.0, 30.0, 50.0}; // [J/kB]

  // Sample spin configurations:
  int n_burnin = 0; // No burn-in yet
  ivec n_cycles = {500, 1000, 8000}; // n monte carlo cycles
  bool Ordered = 0;

  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Open outputfile
  std::ofstream ofile;
  ofile.open("simple_expectation_value_tests.dat");
  ofile << "T   n_cycles   epsilon   m   Cv   X" << std::endl;

  // Compute exp values for different temperatures and n_cycles and write to file
  for (int j=0; j<n_cycles.n_elem; j++)
  {
    for (int i=0; i<T.n_elem; i++)
    {
      double energy;
      double magnetization;
      double heat_capacity;
      double susceptibility;
      expectation_values(L, T(i), Ordered, n_cycles(j), n_burnin, energy, magnetization, heat_capacity, susceptibility);
      ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << n_cycles(j)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << energy
      << std::setw(width) << std::setprecision(decimals) << std::scientific << magnetization
      << std::setw(width) << std::setprecision(decimals) << std::scientific << heat_capacity
      << std::setw(width) << std::setprecision(decimals) << std::scientific << susceptibility
      << std::endl;
    }
  }
  // Close output file
  ofile.close();


  //Problem 5
  L = 20;
  T = {1.0, 2.4}; //[J/kB]
  int n_c = 8000;
  n_burnin = 0; // Still no burn-in yet
  Ordered = 1;

  ofile.open("burnin_test_ordered.dat");
  ofile << "T   n_cycles   epsilon   m" << std::endl;

  // Compute exp values for different numbers of cycles and temperatures and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
     vec en(n_c+1);
     vec mag(n_c+1);
     individual_samples(L, T(i), Ordered, n_c, n_burnin, en, mag); // First element is the random initial state
     // Store the actual samples, not the initial state:
     vec energy = en(span(1, en.n_elem-1));
     vec magnetization = mag(span(1, mag.n_elem-1));
     //cumulative sum of sample energies, divided by number of cycles
     energy = cumsum(energy)/linspace(1, n_c+1, n_c);
     //cumulative sum of sample magnetizations, divided by number of cycles
     magnetization = cumsum(magnetization)/linspace(1, n_c+1, n_c);
     // Add the expectation values after the initial state:
     for (int i=0; i<energy.n_elem; i++)
     {
       en(i+1) = energy(i);
       mag(i+1) = magnetization(i);
     }
     // Write to file :
     for (int j = 0; j < n_c+1; j++){
       ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
       << std::setw(width) << std::setprecision(decimals) << std::scientific << j
       << std::setw(width) << std::setprecision(decimals) << std::scientific << en(j)
       << std::setw(width) << std::setprecision(decimals) << std::scientific << mag(j)
       << std::endl;
     }
   }
  ofile.close();

  // Random start
  Ordered = 0;
  n_c = 8000;
  ofile.open("burnin_test_random.dat");
  ofile << "T   n_cycles   epsilon   m" << std::endl;
  // Compute exp values for different numbers of cycles and temperatures and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
    vec en(n_c+1);
    vec mag(n_c+1);
    individual_samples(L, T(i), Ordered, n_c, n_burnin, en, mag); // First element is the random initial state
    // Store the actual samples, not the initial state:
    vec energy = en(span(1, en.n_elem-1));
    vec magnetization = mag(span(1, mag.n_elem-1));
    //cumulative sum of sample energies, divided by number of cycles
    energy = cumsum(energy)/linspace(1, n_c+1, n_c);
    //cumulative sum of sample magnetizations, divided by number of cycles
    magnetization = cumsum(magnetization)/linspace(1, n_c+1, n_c);
    // Add the expectation values after the initial state:
    for (int i=0; i<energy.n_elem; i++)
    {
      en(i+1) = energy(i);
      mag(i+1) = magnetization(i);
    }
    // Write to file :
    for (int j = 0; j < n_c+1; j++){
      ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << j
      << std::setw(width) << std::setprecision(decimals) << std::scientific << en(j)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << mag(j)
      << std::endl;
   }
  }
  ofile.close();


  // Problem 6
  n_burnin = 3000;
  n_c = 100000; // Number of cycles after burn-in
  Ordered = 0;
  ofile.open("probability_distribution.dat");
  ofile << "T n_cycles epsilon" << std::endl;
  int N=L*L;
  // Get energy samples and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
    vec en(n_c+1);
    vec mag(n_c+1);
    individual_samples(L, T(i), Ordered, n_c, n_burnin, en, mag); // First element is the random initial state
    for (int j = 1; j<n_c+1; j++)
    {
      ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << j
      << std::setw(width) << std::setprecision(decimals) << std::scientific << en(j)
      << std::endl;
    }
  }
  ofile.close();

}
