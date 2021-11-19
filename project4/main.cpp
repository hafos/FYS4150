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
      //std::cout << i <<endl;
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

  //Problem 5
  L = 20;
  T = {1.0, 2.4}; //[J/kB]
  int n_c = 100;
  Ordered = 1;

  ofile.open("burnin_test_ordered.dat");
  ofile << "bounded n_cycles T epsilon m" << std::endl;

  // Compute exp values for different numbers of cycles and temperatures and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
     vec en(n_c);
     vec mag(n_c);
     energy_cycles(L, T(i), Ordered, n_c, n_steps, en);
     magnetization_cycles(L, T(i), Ordered, n_c, n_steps, mag);

     for (int j = 0; j < n_c; j++){
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
  n_c = 200;
  ofile.open("burnin_test_random.dat");
  ofile << "bounded n_cycles T epsilon m" << std::endl;
  // Compute exp values for different numbers of cycles and temperatures and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
    vec en(n_c);
    vec mag(n_c);
    energy_cycles(L, T(i), Ordered, n_c, n_steps, en);
    magnetization_cycles(L, T(i), Ordered, n_c, n_steps, mag);
    for (int j = 0; j < n_c; j++){
      ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << j
      << std::setw(width) << std::setprecision(decimals) << std::scientific << en(j)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << mag(j)
      << std::endl;
   }
  }
  ofile.close();


  //Problem 6
  n_c = 1000;
  ofile.open("probability_distribution.dat");
  ofile << "T n_cycles epsilon" << std::endl;
  // Compute exp values for different numbers of cycles and temperatures and write to file
  for (int i = 0; i < T.n_elem; i++)
  {
    double energy;
    // Do 'n_cycles' monte carlo cycles, store energy after each
    for (int j = 0; j<n_c; j++){
      //Initialize new lattice
      Lattice config(L, T(i), Ordered);
      // Do one monte carlo cycle, length: 'n_steps'
      for (int i = 0; i<n_steps; i++){
        config.spin_flip();
      }
      //get energy
      energy = config.get_energy();
      //write to file
      ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << T(i)
      << std::setw(width) << std::setprecision(decimals) << std::scientific << j
      << std::setw(width) << std::setprecision(decimals) << std::scientific << energy
      << std::endl;
    }
  }
  ofile.close();



}
