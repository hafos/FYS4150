#ifndef __Lattice_hpp__
#define __Lattice_hpp__

#include <armadillo> //vectors and matrices
#include <random>
#include <chrono>

using namespace arma;

class Lattice
{
private:
  //stores the whole lattice
  imat spin_config_;
  //size
  int L_;
  // stores the five possible boltzmann factors
  vec boltzmann_factors_;
  // Mersenne twister random number generator
  std::mt19937 generator_;

public:
  //constructor
  Lattice(int L, double T);

  //return the lattice
  imat get_config();

  //flip one spin
  void spin_flip();

  // Update borders
  void update_borders();

  // Return energy
  double get_energy();

  // Return magnetizaton
  double get_magnetization();
};

#endif
