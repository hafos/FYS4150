#include <iostream>
#include <armadillo>
#include <cmath>

//#include "jacobi.cpp"
//#include "create_tridiag.cpp"

using namespace arma;

void scaling_N(int N, int& iteration){
    // seeing how it scales for eq.6
    double h = 1./(N);
    double a = -1.0/pow(h, 2);
    double d = 2.0/pow(h, 2);
    mat A = create_symmetric_tridiagonal(N, a, d);
    double tol = 1e-6;
    vec val_N(N, fill::zeros);
    mat vec_N(N, N, fill::zeros);
    bool converged;
    jacobi_eigensolver(A, tol, val_N, vec_N,
                       N*10000, iteration, converged);
}

void scaling_dense(int N, int& iteration){
    // seeing how it scales for a random N x N dense matrix
    // generating random symetric matrix of size N x N
    mat A = mat(N, N).randn();
    A = symmatu(A);
    double tol = 1e-6;
    vec val_N(N, fill::zeros);
    mat vec_N(N, N, fill::zeros);
    bool converged;
    jacobi_eigensolver(A, tol, val_N, vec_N,
                       N*1000, iteration, converged);
}
