#include "Schrodinger.hpp"

// Constructor
Schrodinger::Schrodinger(int M_in, double v0, int n_slits)
{
  M_ = M_in;
  // Coordinates:
  x_ = linspace(0,1,M_);
  y_ = linspace<rowvec>(0,1,M_);
  // Set up the potential:
  initialize_potential(v0, n_slits);
}


// Set up an initial state
void Schrodinger::U_init(double xc, double yc, double sx, double sy, double px, double py)
{
  // Center position of wavepacket:     (xc, yc)
  // Standard deviation of wavepacket:  (sx, sy)
  // Momentum of wavepacket:            (px, py)
  // Constructing matrices so that indexing (i,j) corresponds to (x,y)
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
  U_ = exp(exponent);
  // Impose boundary conditions :
  impose_boundaries(U_);
  // Normalize the initial state :
  cx_double n = accu(U_%conj(U_)); // Operator % does element-wise multiplication
  U_ /= sqrt(n);
}

// Impose boundary conditions
void Schrodinger::impose_boundaries(cx_mat& U)
{
  // Set boundaries to 0, only necessary to do this once
  for (int l=0; l<M_; l++)
  {
    U(0, l) = 0;
    U(M_-1, l) = 0;
    U(l, 0) = 0;
    U(l, M_-1) = 0;
  }
}


// Initialize a potential
void Schrodinger::initialize_potential(double v0, int n_slits)
{
  // Parameters of potential :
  // Make it so that these can be read in from a file?
  double wall_thickness = 0.02; // 0.02
  double wall_position = 0.5; // Place at (M_)/2
  double slit_distance = 0.05; // 0.05
  double slit_aperture = 0.05; // 0.05

  // Find start and end of wall in i-index (x):
  int i_0 = index_min(abs( x_ - (wall_position - 0.5*wall_thickness) ));
  int i_1 = index_min(abs( x_ - (wall_position + 0.5*wall_thickness) ));

  // Find starts and ends of wall in j-index (y):
  ivec j(2*n_slits+2);            // To store indices
  j(0)=0;                         // Border
  j(j.n_elem-1)=M_-1;              // Border

  if (n_slits > 0) // If n_slits = 0, j is already filled with the correct indices
  {
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
  }

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
  V_ = V;
}


// Return the wave function
cx_mat Schrodinger::wave_function()
{
  return U_;
}


// Return the potential
sp_mat Schrodinger::potential()
{
  return V_;
}
