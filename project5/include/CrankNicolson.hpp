#ifndef __CrankNicolson_hpp__
#define __CrankNicolson_hpp__

#include <armadillo>
#include <complex>
#include <cmath>


using namespace arma;

// Find the correct vector index from the matrix indexes of u (including the boundaries!)
int vector_index(int M, int i, int j);

// Initialize A and B using M, h, dt and V
void initialize_matrices(int M, double h, double dt, const sp_mat& V,
                          sp_cx_mat& A, sp_cx_mat& B);

#endif
