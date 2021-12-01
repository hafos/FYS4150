#include "CrankNicolson.hpp"


// Find the correct vector index from the matrix indexes of u (including the boundaries!)
int vector_index(int M, int i, int j)
{
  int k = (i-1) + (j-1)*(M-2);
  return k;
}


// Initialize A and B using M, h, dt and V
void initialize_matrices(int M, double h, double dt, const sp_mat& V, sp_cx_mat& A, sp_cx_mat& B)
{
  using namespace std::complex_literals;
  int N = (M-2);
  int N2 = N*N;
  cx_double r = 1i*dt/(2.*h*h);
  cx_vec a(N2);
  cx_vec b(N2);

  // Fill a and b
  for (int k=0; k<N2; k++)
  {
    int j = k/N; // Integer div.
    int i = k - j*N; // from vector_index, and subtract boundary
    a(k) = 1. + 4.*r + 1i*dt/2.*V(i,j);
    b(k) = 1. - 4.*r - 1i*dt/2.*V(i,j);
  }

  // Insert diagonals of a and b
  for (int k=0; k<N2; k++)
  {
    A(k, k) = a(k);
    B(k, k) = b(k);
  }

  // Insert diagonals with only r's
  for (int j=0; j<N2-3; j++)
  {
    A(j, j+3) = -r;
    A(j+3, j) = -r;
    B(j, j+3) = r;
    B(j+3, j) = r;
  }

  // Insert diagonals with r's and some zeros:
  for (int n=0; n<N; n++) // Do each submatrix (there are NxN matrices of NxN elements)
  {
    for (int k=0; k<N-1; k++)
    {
      int kd = k + n*N;
      A(kd, kd+1) = -r;
      A(kd+1, kd) = -r;
      B(kd, kd+1) = r;
      B(kd+1, kd) = r;
    }
  }
}
