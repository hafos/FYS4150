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
  double v0 = 10.0;
  int n_slits = 2;
  Schrodinger syst1(M, v0, n_slits);
  sp_mat V = syst1.potential(); // The potential in Schrodinger class is (MxM)
  sp_cx_mat A(N2, N2);
  sp_cx_mat B(N2, N2);

  initialize_matrices(M, h, dt, V, A, B);

  std::cout << "A real + imag :" << endl << real(cx_mat(A)) << endl << imag(cx_mat(A)) << endl;
  std::cout << "B real + imag :" << endl << real(cx_mat(B)) << endl << imag(cx_mat(B)) << endl;


  std::cout << "Test of solver :" << endl;
  cx_vec u = {1,0,0,  0,0,0,  0,0,0};
  cx_vec u_new = compute_next_step(u, A, B);
  std::cout << u_new << endl;


  M = 5;
  N = M-2;
  N2 = N*N;
  std::cout << "Test of system initialization: " << endl;
  Schrodinger syst(M, v0, n_slits);
  syst.U_init(0.5, 0.5, 0.3, 0.3, 0.1, 0.1);
  std::cout << real(syst.wave_function()%conj(syst.wave_function())) << endl;


  std::cout << "Test of potential initialization: " << endl;
  //syst.initialize_potential(); // Is done in constructor...
  std::cout << mat(syst.potential()) << endl;
  // If the slits are not visible here it is because;
  // - M is too small
  // - slit_aperture is too small for the M
  // Will of course be better when simulating with small h -> large M.
}
