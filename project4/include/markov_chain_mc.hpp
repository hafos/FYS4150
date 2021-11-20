#ifndef __markov_chain_mc_hpp__
#define __markov_chain_mc_hpp__

#include <iostream>
#include <armadillo> // Vectors and matrices
#include <cmath>

#include "Lattice.hpp"

using namespace arma;

// Compute expectation values
void expectation_values(int L, double T, bool Ordered, int n_cycles,
                        double& eps, double& mag, double& heat, double& sus);

// Compute energy and magnetization expectation value for each cycle
void energy_cycles(int L, double T, bool Ordered, int n_cycles, vec&e, vec&m);

#endif
