#include <iostream>
#include <armadillo>
#include <cmath>

#include "ilode.cpp"

using namespace arma;

void jacobi_rotate(arma::mat& B, arma::mat& R, int k, int l){
  // Make idetity matrix S:
  mat S(B.n_cols, B.n_cols, fill::zeros);
  for (int i=0; i<B.n_cols; i++){
    S(i,i) = 1;
  }

  // Get k, l, tau, t, c, s:
  double maxval = max_offdiag_symmetric(B, k, l); // sets new k and l
  double tau = (B(l,l) - B(k,k))/(2*B(k,l));
  vec ts(2);
  ts(0) = - tau + sqrt(1 + pow(tau, 2));
  ts(1) = - tau - sqrt(1 + pow(tau, 2));
  double t = min(ts);
  double c = 1/sqrt(1+pow(t,2));
  double s = t*c;

  // Update S:
  S(k, k) = c;
  S(l, l) = c;
  S(k, l) = -s;
  S(l, k) = s;

  // Update B:
  B = (S.t()*B)*S;
  // Update R:
  R = R*S;
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
    jacobi_rotate(B, R, k, l); // This updates B and R
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
    for (int j=0; j<A.n_cols; j++){
      eigenvectors(i,j) = R(i,j);
    }
  }
}
