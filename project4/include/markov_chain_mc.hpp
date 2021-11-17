#ifndef __markov_chain_mc_hpp__
#define __markov_chain_mc_hpp__

#include <iostream>
#include <armadillo> // Vectors and matrices
#include <cmath>

#include "Lattice.hpp"

using namespace arma;

// Compute expectation values
void expectation_values(int L, double T, bool Ordered, int n_cycles, int n_steps,
                        double& eps, double& mag, double& heat, double& sus);


#endif
