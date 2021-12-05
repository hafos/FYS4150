#include "Schrodinger.hpp"

// Constructor
Schrodinger::Schrodinger(int M_in)
{
  M_ = M_in;
  x_ = linspace(0,1,M_);
  std::cout << x_ << endl;
  y_ = linspace<rowvec>(0,1,M_);
  V_ = initialize_potential();
}


// Set up an initial state
cx_mat Schrodinger::u_init(double xc, double yc, double sx, double sy, double px, double py)
{
  // Constructing matrices so that indexing (i,j) corresponds to (x,y)
  // I have double checked this by printing elements X and Y :)
  // Printing the entire matrices puts i in vertical (top to bottom) and j in horizontal (left to right)
  // Make coordinate grid :
  mat X(M_, M_);
  mat Y(M_, M_);
  X.each_col() = x_;
  Y.each_row() = y_;
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
  for (int l=0; l<M_; l++) // This works fine..
  {
    u_init(0, l) = 0;
    u_init(M_-1, l) = 0;
    u_init(l, 0) = 0;
    u_init(l, M_-1) = 0;
  }
  // The following might be a more elegant way of doing it but it doesn't work..
  //using namespace std::complex_literals;
  //uvec indices = {0, M_-1};
  //u_init.each_col(indices) = 0i;
  //u_init.each_row(indices) = 0i;
}


// Initialize a potential
sp_mat Schrodinger::initialize_potential()
{
  // Parameters of potential :
  // Make it so that these can be read in from a file?
  double wall_thickness = 0.02; // 0.02
  double wall_position = 0.5; // Place at (M_)/2
  double slit_distance = 0.05; // 0.05
  double slit_aperture = 0.05; // 0.05
  double v0 = 10; // Strength of potential in wall..
  int n_slits = 2; // Number of slits

  // Find start and end of wall in i-index (x):
  int i_0 = index_min(abs( x_ - (wall_position - 0.5*wall_thickness) ));
  int i_1 = index_min(abs( x_ - (wall_position + 0.5*wall_thickness) ));

  // Find starts and ends of wall in j-index (y):
  ivec j(2*n_slits+2);            // To store indices
  j(0)=0;                         // Border
  j(j.n_elem-1)=M_-1;              // Border

  // Set a few characteristic lengths
  double midpoint = 0.5*(max(y_) - min(y_)); // or y_(-1) - y_(0)
  double length = n_slits*slit_aperture + (n_slits - 1)*slit_distance; // From start of first slit to the end of the last
  double position = (midpoint - length/2); // Start from start of first slit

  // Set the first internal wall border :
  j(1) = index_min(abs( y_ - position ) );

  // Now loop over internal wall borders after the first one :
  for (int l=2; l<j.n_elem-2; l++)
  {
    position += slit_aperture;
    j(l) = index_min(abs( y_ - position ) );
    position += slit_distance;
    j(l+1) = index_min(abs( y_ - position ) );
    l += 1; // Skip one
  }

  // Set the last internal wall border :
  j(j.n_elem-2) = index_min(abs( y_ - (position + slit_aperture) ) );

  // Rows: from i0 to i1 is wall
  // Columns: From j0 to j1 and from j2 to j3 and from j4 to j5 etc is wall.

  sp_mat V(M_, M_);

  // Loop over the indices and set wall areas in potential :
  for (int l=0; l<j.n_elem; l++)
  {
    mat sec(i_1-i_0+1, j(l+1)-j(l)+1);
    sec.fill(v0); // Doing this as Ubuntu has arma v.9.8 and not >10.6 pre-installed -> fill::value(v0) not available
    V( span(i_0, i_1), span( j(l), j(l+1) ) ) = sec;
    l += 1; // Skip one
  }
  return V;
}
