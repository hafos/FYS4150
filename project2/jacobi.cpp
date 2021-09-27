#include <iostream>
#include <armadillo>
#include <cmath>

#include "ilode.cpp"

using namespace arma;

void jacobi_rotate(arma::mat& B, arma::mat& R){
  int k = 0;
  int l = 0;

  // Get k, l, tau, t, c, s:
  double maxval = max_offdiag_symmetric(B, k, l); // sets new k and l
  if(k == l){
    std::cout << "Something went wrong: k = l" << endl;
  }
  double tau = (B(l, l) - B(k, k))/(2*B(k, l));
  double t;
  if(tau > 0){
    t = - tau + sqrt(1 + pow(tau, 2));
  }
  if(tau <= 0){
    t = - tau - sqrt(1 + pow(tau, 2));
  }
  double c = 1/sqrt(1 + pow(t,2));
  double s = t*c;

  // Update B:
  //B = (S.t()*B)*S;
  double B_kk = B(k, k)*pow(c, 2) - 2*B(k, l)*c*s + B(l, l)*pow(s, 2);
  B(l, l) = B(l, l)*pow(c, 2) + 2*B(k, l)*c*s + B(k, k)*pow(s, 2);
  B(k, k) = B_kk;
  B(k, l) = 0;
  B(l, k) = 0;
  double B_ik = 0;
  for (int i=0; i < B.n_cols; i++){
    if(i != k){
      if(i != l){
        B_ik = B(i, k)*c - B(i, l)*s;
        B(k, i) = B_ik;
        B(i, l) = B(i, l)*c + B(i, k)*s;
        B(l, i) = B(i, l);
        B(i, k) = B_ik;
      }
    }
  }
  // Update R:
  //R = R*S;
  double R_ik;
  for (int i=0; i < B.n_cols; i++){
    R_ik = R(i, k)*c - R(i, l)*s;
    R(i, l) = R(i, l)*c + R(i, k)*s;
    R(i, k) = R_ik;
  }
}

void jacobi_eigensolver(const arma::mat& A, double tolerance, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged){
  iterations = 0;
  converged = 0;

  // R = I
  mat R(A.n_cols, A.n_cols, fill::zeros);
  for (int i=0; i < A.n_cols; i++){
    R(i, i) = 1;
  }
  // B = A
  mat B(A.n_cols, A.n_cols, fill::zeros);
  for (int i=0; i < A.n_cols; i++){
    for (int j=0; j < A.n_cols; j++){
      B(i, j) = A(i, j);
    }
  }

  int k=0;
  int l=0;
  double maxoff = max_offdiag_symmetric(B, k, l);
  while(pow(maxoff, 2) > tolerance){
    jacobi_rotate(B, R); // This updates k, l then B and R
    iterations = iterations + 1;
    double maxoff = max_offdiag_symmetric(B, k, l);
    // If-tests for convergence etc.:
    if(iterations >= maxiter){
      std::cout << "Jacobi rotation method not converged after " << iterations;
      std::cout << " iterations. (tolerance = " << tolerance << " )" << endl;
      break;
    }
    if(pow(maxoff, 2) <= tolerance){
      converged = 1;
      break;
    }
  }
  for (int i=0; i < A.n_cols; i++){
    eigenvalues(i) = B(i,i);
    for (int j=0; j < A.n_cols; j++){
      eigenvectors(i, j) = R(i, j);
    }
  }
}
