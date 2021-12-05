#include "Schrodinger.hpp"

// Constructor
Schrodinger::Schrodinger(int M_in)
{
  M = M_in;
  x = linspace(0,1,M);
  std::cout << x << endl;
  y = linspace<rowvec>(0,1,M);
}


// Set up an initial state
cx_mat Schrodinger::u_init(double xc, double yc, double sx, double sy, double px, double py)
{
  // Constructing matrices so that indexing (i,j) corresponds to (x,y)
  // I have double checked this by printing elements X and Y :)
  // Printing the entire matrices puts i in vertical (top to bottom) and j in horizontal (left to right)
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


// Initialize a potential
sp_mat Schrodinger::initialize_potential()
{
  // Make it so that these can be read in from a file?
  double wall_thickness = 0.02; // 0.02
  double wall_position = 0.5; // Place at (M)/2
  double slit_distance = 0.05; // 0.05
  double slit_aperture = 0.05; // 0.05
  double v0 = 10; // Strength of potential in wall..

  // Find start and end of wall in i-index (x):
  int i_0 = index_min(abs( x - (wall_position - 0.5*wall_thickness) ));
  int i_1 = index_min(abs( x - (wall_position + 0.5*wall_thickness) ));

  // Find starts and ends of wall in j-index (y):
  int j_0 = 0;
  int j_1 = index_min(abs( y - ( 0.5*(max(y) - min(y)) - 0.5*slit_distance - slit_aperture ) ));
  int j_2 = index_min(abs( y - ( 0.5*(max(y) - min(y)) - 0.5*slit_distance ) ));
  int j_3 = index_min(abs( y - ( 0.5*(max(y) - min(y)) + 0.5*slit_distance ) ));
  int j_4 = index_min(abs( y - ( 0.5*(max(y) - min(y)) + 0.5*slit_distance + slit_aperture ) ));
  int j_5 = M-1;

  // Rows: from i0 to i1
  // Columns: From j0 to j1 and from j2 to j3 and from j4 to j5.

  sp_mat V(M, M);
  mat sec1(i_1-i_0+1, j_1-j_0+1);
  sec1.fill(v0); // Doing this as Ubuntu has arma v.9.8 and not >10.6 pre-installed -> fill::value(v0) not available
  mat sec2(i_1-i_0+1, j_3-j_2+1);
  sec2.fill(v0);
  mat sec3(i_1-i_0+1, j_5-j_4+1);
  sec3.fill(v0);
  V( span(i_0, i_1), span(j_0, j_1) ) = sec1;
  V( span(i_0, i_1), span(j_2, j_3) ) = sec2;
  V( span(i_0, i_1), span(j_4, j_5) ) = sec3;
  return V;
}
