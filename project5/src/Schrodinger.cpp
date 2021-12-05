#include "Schrodinger.hpp"

// Constructor
Schrodinger::Schrodinger(int M_in)
{
  M = M_in;
  x = linspace(0,1,M);
  //std::cout << x << endl;
  y = linspace<rowvec>(0,1,M);
}


// Set up an initial state
cx_mat Schrodinger::u_init(double xc, double yc, double sx, double sy, double px, double py)
{
  // Constructing matrices so that indexing (i,j) corresponds to (x,y)
  // I have double checked this by printing elements X and Y :)
  // Printing the entire matrices puts i in vertical and j in horizontal
  // Make coordinate grid :
  mat X(M, M);
  mat Y(M, M);
  X.each_col() = x;
  Y.each_row() = y;
  // Make grid of real and imaginary part of exponent :
  mat rel = - pow( (X - xc), 2)/(2.*sx*sx) - pow( (Y - yc), 2)/( 2.*sy*sy );
  mat imag = px*(X - xc) + py*(Y - yc);
  // Put the exponent together as a complex matrix :
  cx_mat exponent(rel, imag);
  // Compute the initial state grid :
  cx_mat u_in = exp(exponent);
  // Impose boundary conditions :
  impose_boundaries(u_in);
  // Normalize the initial state :
  cx_double n = accu(u_in%conj(u_in)); // Operator % does element-wise multiplication
  u_in /= sqrt(n);
  n = accu(u_in%conj(u_in)); // Operator % does element-wise multiplication
  //std::cout << "Normalized if 1 : " << n << endl;
  return u_in;
}

// Impose boundary conditions on initial state
void Schrodinger::impose_boundaries(cx_mat& u_init)
{
  for (int l=0; l<M; l++) // This works fine..
  {
    u_init(0, l) = 0;
    u_init(M-1, l) = 0;
    u_init(l, 0) = 0;
    u_init(l, M-1) = 0;
  }
  // The following might be a more elegant way of doing it but it doesn't work..
  //using namespace std::complex_literals;
  //uvec indices = {0, M-1};
  //u_init.each_col(indices) = 0i;
  //u_init.each_row(indices) = 0i;
}
