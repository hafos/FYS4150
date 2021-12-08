#include "CrankNicolson.hpp"


// Find the correct vector index from the matrix indexes of u (including the boundaries!)
int vector_index(int M, int i, int j)
{
  int k = (i-1) + (j-1)*(M-2);
  return k;
}

// Find the correct matrix index from the vector indexes of u (including the boundaries!)
void matrix_index(int M, int k, int& i, int& j)
{
  // Both methods give same answer ;
  //i = k % (M-2) + 1; //% is the modulus here
  //j = (k-i+1)/(M-2) + 1;
  j = k/(M-2) + 1; // Integer div., for V with boundaries, size (MxM)
  i = (k + 1) - (j - 1)*(M-2); // from vector_index
}


// Initialize A and B using M, h, dt and V
void initialize_matrices(int M, double h, double dt, const sp_mat& V, sp_cx_mat& A, sp_cx_mat& B)
{
  // V of size (MxM)
  // A and B of size ((M-2)**2 x (M-2)**2)
  using namespace std::complex_literals; // For 1i notation
  int N = (M-2);
  int N2 = N*N;
  cx_double r = 1i*dt/(2.*h*h);
  cx_vec a(N2);
  cx_vec b(N2);

  // Fill a and b
  for (int k=0; k<N2; k++)
  {
    int i, j;
    matrix_index(M, k, i, j);
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

// Compute and return the next step u^(n+1)
cx_vec compute_next_step(const cx_vec& u, const sp_cx_mat& A, const sp_cx_mat& B)
{
  // Compute the right hand side :
  cx_vec b = B*u;

  // Now solve the matrix equation Au^(n+1) = b
  cx_vec u_new = spsolve(A, b, "superlu");

  return u_new;
}

// Convert matrix to vector representation
cx_vec mat2vec(int M, cx_mat U)
{
  int N = M-2;
  int N2 = N*N;

  cx_vec u_vec(N2); //initialize

  for (int j = 1; j < M-1; j++)
  {
    for (int i = 1; i < M-1; i++)
    {
      int k = vector_index(M, i, j);
      //std::cout << "i: " << i << "  j: " << j << "  -->  k: " << k << endl;
      u_vec(k) = U(i, j);
    }
  }
  return u_vec;
}

// Convert vector to matrix representation
cx_mat vec2mat(int M, cx_vec u)
{
  cx_mat U_mat(M, M, fill::zeros); //initialize
  for (int k = 0; k < u.n_rows; k++)
  {
    int i, j;
    matrix_index(M, k, i, j);
    //std::cout << "k = " << k << " --> i: " << i << "  j: " << j << endl;
    U_mat(i, j) = u(k);
  }
  return U_mat;
}
