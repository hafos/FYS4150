#include <iostream>
#include <armadillo>
#include <cmath>

#include "ilode.cpp"

using namespace arma;

void jacobi_rotate(arma::mat& B, arma::mat& eigenvectors, arma::mat& R, int k, int l){
  // Reset R to identity matrix
  R(k, l) = 0;
  R(l, k) = 0;
  R(k, k) = 1;
  R(l, l) = 1;

  // Get k, l, tau, t, c, s:
  double maxval = max_offdiag_symmetric(B, k, l); // THIS SETS NEW k AND l
  double tau = (B(l,l) - B(k,k))/(2*B(k,l));
  vec ts(2);
  ts(0) = - tau + sqrt(1 + pow(tau, 2));
  ts(1) = - tau - sqrt(1 + pow(tau, 2));
  double t = min(ts);
  double c = 1/sqrt(1+pow(t,2));
  double s = t*c;

  // Update R:
  R(k, k) = c;
  R(l, l) = c;
  R(k, l) = -s;
  R(l, k) = s;

  // Update B:
  B = (R.t()*B)*R;
  // Update eigenvectors:
  eigenvectors = eigenvectors*R;
}

void jacobi_eigensolver(const arma::mat& A, double tolerance, arma::vec& eigenvalues, arma::mat& eigenvectors,
                        const int maxiter, int& iterations, bool& converged){
  //double tolerance = 1e-8;
  iterations = 0;
  converged = 0;

  int k = 0;
  int l = 1;
  // Make R and eigenvectors identity matrices to edit later
  mat R(A.n_cols, A.n_cols, fill::zeros);
  for (int i=0; i<A.n_cols; i++){
    eigenvectors(i,i) = 1;
    R(i,i) = 1;
  }
  // Copy A to B
  mat B(A.n_cols, A.n_cols, fill::zeros);
  for (int i=0; i<A.n_cols; i++){
    for (int j=0; j<A.n_cols; j++){
      B(i,j) = A(i,j);
    }
  }

  while(pow(B, 2).max() > tolerance){
    jacobi_rotate(B, eigenvectors, R, k, l);
    iterations = iterations + 1;

    // If-tests for convergence etc.:
    if(iterations > maxiter){
      std::cout << "Not converged " << endl;
      break;
    }
    if(pow(B, 2).max() <= tolerance){
      std::cout << "Converged! " << endl;
      converged = 1;
    }
    if(k==l){
      std::cout << "Something went wrong: k = l" << endl;
      break;
    }
  }
  std::cout << "Loop Finished " << endl;
  for (int i=0; i<A.n_cols; i++){
    eigenvalues(i) = B(i,i);
  }
}
