#include <iostream>
#include <iomanip> //Writing to file
#include <armadillo> //Vectors and matrices
#include <cmath>

#include "Lattice.hpp"

using namespace arma;

int main()
{
  //some test code for the Lattice methods
  int L = 2;
  Lattice config(L);
  std::cout << config.get_config() << endl;

  config.spin_flip();
  std::cout <<config.get_config() << endl;
}
