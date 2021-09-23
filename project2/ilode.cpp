#include <iostream>
#include <armadillo>
#include <cmath>

using namespace arma;

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l){
  double max = 0;
  for (int i = 1; i < A.n_rows; i++){
    for (int j = 0; j < i; j++){
      if (A(i,j) > max){
        max = A(i,j);
        k = i;
        l = j;
      }
    }
  }
  return max;
}

int main() {
  // Make test matrix:
  mat A(4, 4, fill::zeros);
  for (int i = 0; i < 4; i++){
      A(i, i) = 1;
      if(i==0 || i==3 ){A(i, 3-i) = 0.5;
      } else{A(i, 3-i) = -0.7 ;
      }
    }

  // Test the function:
  int k = 0;
  int l = 0;
  double m = max_offdiag_symmetric(A, k, l);

  std::cout << endl << "Test of max_offdiag_symmetric:" << endl;
  std::cout << "k: "<< k << "  l: " << l << "  max: " << m << endl;
  return 0;
}
