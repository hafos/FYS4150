#include <iostream> // Printing to terminal and file

#include "CrankNicolson.hpp"
#include "Schrodinger.hpp"

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

  std::cout << "A :" << endl << cx_mat(A) << endl  << "B :" << endl << cx_mat(B) << endl;


  std::cout << "Test of solver :" << endl;
  cx_vec u = {1,0,0,  0,0,0,  0,0,0};
  cx_vec u_new = compute_next_step(u, A, B);
  std::cout << u_new << endl;


  M = 10;
  std::cout << "Test of system initialization: " << endl;
  Schrodinger syst(M);
  cx_mat uin;
  uin = syst.u_init(1., 1., 1., 1., 1., 1.);
  std::cout << uin << endl;


  std::cout << "Test of potential initializer: " << endl;
  sp_mat V2;
  V2 = syst.initialize_potential();
  std::cout << mat(V2) << endl;
  // If the slits are not visible here it is because;
  // - M is too small
  // - slit_aperture is too small
}
