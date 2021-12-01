#include <iostream> // Printing to terminal and file

#include "CrankNicolson.hpp"

int main()
{
  std::cout << "Test of vector_index :" << endl;
  int M = 5;
  ivec i = {1,2,3};
  ivec j = {1,2,3};
  for (int ij=0; ij<3; ij++)
  {
    for(int ii=0; ii<3; ii++)
    {
      int k = vector_index(M, i(ii), j(ij));
      std::cout << "i: " << i(ii) << "  j: " << j(ij) << "  -->  k: " << k << endl;
    }
  }
  std::cout << endl;


  std::cout << "Test of initialize_matrices :" << endl;
  double h = 1;
  double dt = 1;
  int N = M-2;
  int N2 = N*N;
  sp_mat V(N,N);
  sp_cx_mat A(N2, N2);
  sp_cx_mat B(N2, N2);

  initialize_matrices(M, h, dt, V, A, B);

  std::cout << "A :" << endl << A << endl  << "B :" << endl << B << endl;
}
