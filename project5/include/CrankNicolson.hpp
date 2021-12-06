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

// Compute and return the next step u^(n+1)
cx_vec compute_next_step(const cx_vec& u, const sp_cx_mat& A, const sp_cx_mat& B);

// Convert matrix to vector representation
cx_vec mat2vec(int M, cx_mat U);

// Convert vector to matrix representation
cx_mat vec2mat(int M, cx_vec u);

#endif
