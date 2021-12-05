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
  vec x;
  rowvec y;

  // Number of points:
  int M;

public:
  // Constructor
  Schrodinger(int M);

  // Set up an initial state
  cx_mat u_init(double xc, double yc, double sx, double sy,
                double px, double py);
                
  // Impose boundary conditions on initial state
  void impose_boundaries(cx_mat& u_init);
};



#endif
