#include <iomanip> // Writing to file
#include <string> // Turn to string
#include <sstream> // String formatting

#include "Lattice.hpp"
#include "timecheck.hpp"
#include "markov_chain_mc.hpp"

using namespace arma;

int main()
{
  // Problem 4:
  int L = 2;
  vec T = {0.5, 1.0, 1.5, 2.5, 4.5, 7.0, 10.0, 30.0, 50.0}; // [J/kB]
  int n_chains = 10; // Run several chains

  // Sample spin configurations:
  ivec n_cycles = {2000, 5000, 8000}; // n monte carlo cycles
  bool Ordered = 0;

  // spacing in outputfile
  int width = 21; int decimals = 9;

  // Open outputfile
  std::ofstream ofile;
  ofile.open("simple_expectation_value_tests.dat");
  ofile << "T   n_cycles   epsilon   m   Cv   X" << std::endl;

  // Compute exp values for different temperatures and n_cycles/steps and write to file
  for (int j=0; j<n_cycles.n_elem; j++)
  {
    for (int i=0; i<T.n_elem; i++)
    {
      double avg_e = 0;
      double avg_m = 0;
      double avg_C = 0;
      double avg_X = 0;
      for (int k=0; k<n_chains; k++)
      {
        double energy;
        double magnetization;
        double heat_capacity;
        double susceptibility;
        expectation_values(L, T(i), Ordered, n_cycles(j), energy, magnetization, heat_capacity, susceptibility);
        avg_e += energy;
        avg_m += magnetization;
        avg_C += heat_capacity;
        avg_X += susceptibility;
      }
      ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << n_cycles(j)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << avg_e/n_chains
      << std::setw(width) << std::setprecision(decimals) << std::scientific << avg_m/n_chains
      << std::setw(width) << std::setprecision(decimals) << std::scientific << avg_C/n_chains
      << std::setw(width) << std::setprecision(decimals) << std::scientific << avg_X/n_chains
      << std::endl;
    }
  }
  // Close output file
  ofile.close();

  // Testing how fast methods are
  //timecheck();


  //Problem 5
  L = 20;
  T = {1.0, 2.4}; //[J/kB]
  int n_c = 8000;
  Ordered = 1;

  ofile.open("burnin_test_ordered.dat");
  ofile << "T   n_cycles   epsilon   m" << std::endl;

  // Compute exp values for different numbers of cycles and temperatures and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
     vec en(n_c+1);
     vec mag(n_c+1);
     energy_cycles(L, T(i), Ordered, n_c, en, mag);
     for (int j = 0; j < n_c+1; j++){
       ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
       << std::setw(width) << std::setprecision(decimals) << std::scientific << j
       << std::setw(width) << std::setprecision(decimals) << std::scientific << en(j)
       << std::setw(width) << std::setprecision(decimals) << std::scientific << mag(j)
       << std::endl;
     }
   }
  ofile.close();

  //random start
  Ordered = 0;
  n_c = 8000;
  ofile.open("burnin_test_random.dat");
  ofile << "T   n_cycles   epsilon   m" << std::endl;
  // Compute exp values for different numbers of cycles and temperatures and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
    vec en(n_c+1);
    vec mag(n_c+1);
    energy_cycles(L, T(i), Ordered, n_c, en, mag);
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
  n_c = 3000; // Number of cycles after burn-in
  Ordered = 0;
  ofile.open("probability_distribution.dat");
  ofile << "T n_cycles epsilon" << std::endl;
  int N=L*L;
  // Compute exp values for different numbers of cycles and temperatures and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
   double energy;
   // Do 'n_cycles' monte carlo cycles, store energy after each
   // Initialize new lattice
   Lattice config(L, T(i), Ordered);
   // Burn-in:
   for (int j = 0; j<1000; j++){
     for (int k = 0; k<N; k++){
       config.spin_flip();
     }
   }
   for (int j = 0; j<n_c; j++){
     // Do one monte carlo cycle, length: N
     for (int k = 0; k<N; k++){
       config.spin_flip();
     }
     //get energy
     energy = config.get_energy();
     //write to file
     ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
     << std::setw(width) << std::setprecision(decimals) << std::scientific << j+1
     << std::setw(width) << std::setprecision(decimals) << std::scientific << energy/N
     << std::endl;
   }
  }
  ofile.close();



}
