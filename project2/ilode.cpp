#include <iostream>
#include <armadillo>
#include <cmath>

using namespace arma;

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
  double max = 0;
  for (int i = 1; i < A.n_cols; i++){
    for (int j = 0; j < i; j++){
      if (std::abs(A(i,j)) > max){
        max = std::abs(A(i,j));
        k = i;
        l = j;
      }
    }
  }
  return max;
}

void test_ilode() {
  // Make test matrix:
  mat test_matrix(4, 4, fill::zeros);
  for (int i = 0; i < 4; i++){
      test_matrix(i, i) = 1;
      if(i==0 || i==3 ){test_matrix(i, 3-i) = 0.5;
      } else{test_matrix(i, 3-i) = -0.7 ;
      }
    }

  // Test the function:
  int k = 0;
  int l = 0;
  double m = max_offdiag_symmetric(test_matrix, k, l);

  std::cout << "Test of max_offdiag_symmetric() :" << endl;
  std::cout << " k: "<< k << "  l: " << l << "  max: " << m << endl;
  std::cout << " The  test matrix: " << endl;
  test_matrix.print();
  std::cout << endl;
}
