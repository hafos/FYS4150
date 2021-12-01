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



  // double h = 1;
  // double dt = 1;
  // sp_mat V;
  // sp_cx_mat A;
  // sp_cx_mat B;

  // initialize_matrices(M, h, dt, V, A, B);

  // std::cout << A << endl << endl << B << endl;
}
