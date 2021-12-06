#ifndef __Schrodinger_hpp__
#define __Schrodinger_hpp__

#include <armadillo>
#include <complex>
#include <cmath>

using namespace arma;

class Schrodinger
{
private:
  // Grid
  vec x_;
  rowvec y_;

  // Number of points
  int M_;

  // Potential
  sp_mat V_;

  // The wave function:
  cx_mat U_;

public:
  // Constructor
  Schrodinger(int M_in, double v0, int n_slits);

  // Set up an initial state
  void U_init(double xc, double yc, double sx, double sy,
                double px, double py);

  // Impose boundary conditions
  void impose_boundaries(cx_mat& U);

  // Initialize a potential
  void initialize_potential(double v0, int n_slits);

  // Return the wave function
  cx_mat wave_function();

  // Return the potential
  sp_mat potential();
};



#endif
