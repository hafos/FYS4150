#include <iostream>
#include <iomanip> //Writing to file
#include <armadillo> //Vectors and matrices
#include <cmath>

#include "Lattice.hpp"

using namespace arma;

void expval_MCMC(int L, double T, int n_cycles, int n_steps, double& eps, double& mag, double& heat, double& sus){
  double energy;
  double E2;
  double E;
  double magnetization;
  double M2;
  double M;

  int N = L*L;

  for (int j=0; j<n_cycles; j++){
    Lattice config(L, T);
    for (int i=0; i<n_steps; i++){
      config.spin_flip();
    }
    energy = config.get_energy();
    magnetization = abs(config.get_magnetization());
    E += energy/n_cycles;
    M += magnetization/n_cycles;
    E2 += energy*energy/n_cycles;
    M2 += magnetization*magnetization/n_cycles;
  }
  // Calculate the values:
  eps = E/N; // Units [J]
  mag = M/N; // Units [1] or [spin]
  heat = (E2 - E*E)/(T*T*N); // units [kB]
  sus = (M2 - M*M)/(T*N); // units [1/J] or [spin^2 / J]
}

int main()
{
  //some test code for the Lattice methods
  int L = 2;
  double T = 1; // [J/kB]
  //Lattice config(L, T);
  //std::cout << config.get_config() << endl;

  //config.spin_flip();
  //std::cout <<config.get_config() << endl;

  //config.spin_flip();
  //std::cout <<config.get_config() << endl;

  // Sample spin configurations:
  int n_cycles = 200; // n monte carlo cycles
  int n_steps = 50; // steps in each monte carlo cycle

  double energy;
  double magnetization;
  double heat_capacity;
  double susceptibility;

  expval_MCMC(L, T, n_cycles, n_steps, energy, magnetization, heat_capacity, susceptibility);

  std::cout << " <eps> " << energy << "  J" << endl;
  std::cout << " <|m|> " << magnetization << endl;
  std::cout << " Cv " << heat_capacity << "  kB" << endl;
  std::cout << " X " << susceptibility << "  1/J" << endl;
}
