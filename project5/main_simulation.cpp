#include <iostream> // Printing to terminal and file

#include "CrankNicolson.hpp"
#include "Schrodinger.hpp"

using namespace arma;

int main(int argc, char* argv[])
{
  // Check number of command-line arguments

  if (argc != 11)  // Expect 10 command-line arguments
  {
    // Get the name of the executable file
    std::string executable_name = argv[0];

    std::cerr << "Error: Wrong number of input arguments." << std::endl;
    std::cerr << "Usage: " << executable_name << " <some integer> <some double> <some string>" << std::endl;

    // Exit program with non-zero return code to indicate a problem
    return 1;
  }

  // read in parameters from command line
  double h_in = atof(argv[1]); // position step size (not Planck constant)
  double dt_in = atof(argv[2]); // time step size
  double Tmax_in = atof(argv[3]); // end of time interval
  double xc_in = atof(argv[4]); // center of initial wave packet (x)
  double sx_in = atof(argv[5]); // initial width of wave packet (x)
  double px_in = atof(argv[6]); // initial momentum of wave packet (x)
  double yc_in = atof(argv[7]); // center of initial wave packet (y)
  double sy_in = atof(argv[8]); // initial width of wave packet (y)
  double py_in = atof(argv[9]); // initial momentum of wave packet (y)
  double v0_in = atof(argv[10]); // potential barrier height

  int M = 5; // just for testing
  int N = M-2;
  int N2 = N*N;
  int n_slits = 2;
  Schrodinger syst(M, v0_in, n_slits); // sets up potential and system

  syst.U_init(xc_in, yc_in, sx_in, sy_in, px_in, py_in); // sets up initial state matrix
  cx_mat U = syst.wave_function(); // matrix form
  //std::cout << U << endl;
  cx_vec u_vec = mat2vec(M, U); // vector form

  sp_mat V = syst.potential();
  sp_cx_mat A(N2, N2);
  sp_cx_mat B(N2, N2);
  initialize_matrices(M, h_in, dt_in, V, A, B); // sets up Crank-Nicolson

  // loop over time
  // store initial condition as first slice of probability density cube
  cube P_n(M, M, 1, fill::zeros);
  P_n.slice(0) = real(U%conj(U));

  double t = 0;
  mat p(M, M);
  while (t < Tmax_in)
  {
    u_vec = compute_next_step(u_vec, A, B);
    U = vec2mat(M, u_vec);
    p = real(U%conj(U));
    P_n = join_slices(P_n, p);
    t = t+dt_in;
  }

  std::cout << P_n.tail_slices(1) << endl;

  // save as binary data file
  P_n.save("schrodinger_probability.bin");


}
