#ifndef __Lattice_hpp__
#define __Lattice_hpp__

#include <armadillo> //vectors and matrices
#include <random>
#include <chrono>

using namespace arma;

class Lattice
{
private:
  // Stores the whole lattice
  imat spin_config_;

  // Size
  int L_;

  // Stores the five possible boltzmann factors
  vec boltzmann_factors_;

  // Stores the five possible energy differences
  vec delta_E_;

  // Mersenne twister random number generator
  std::mt19937 generator_;

public:
  // Constructor
  Lattice(int L, double T, bool Ordered);

  // Return the lattice
  imat get_config();

  // Pick one random spin. Flip it if accepted
  void spin_flip();

  // Update all borders
  void update_borders();

  // Return total energy
  double get_energy();

  // Return total magnetizaton
  double get_magnetization();

  // Stores the current total energy
  double total_E_;

  // Stores the current total magnetization
  double total_M_;
};

#endif
