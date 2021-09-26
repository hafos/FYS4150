//
//  On Mac, compile in the terminal using command:
//  g++ main.cpp -std=c++11 -o main.exe -larmadillo
//

#include <iostream>
#include <armadillo>
#include <cmath>
#include "jacobi.cpp" // includes ilode.cpp
#include "vector_check.cpp"

using namespace arma;


int main() {
    //Problem 3: initialize a tridiagonal 6x6 matrix A
    //with diagonal 2/h^2 and sub/superdiagonal -1/h^2
    //h = 1/(N+1) so here it's 1/7
    int N = 6;
    mat A(6, 6, fill::zeros);
    double a = -1.0/pow(1./7., 2);
    double d = 2.0/pow(1./7., 2);

    for (int i = 0; i < 6; i++){
        A(i, i) = d;
        if (i+1 < 6){
            A(i, i+1) = a;
            A(i+1, i) = a;
        }
    }

    // Solve the eigenvalue problem numerically:
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);
    eigvec = normalise(eigvec);

    // Compute analytical result:
    vec val_ana(6);
    mat vec_ana(6, 6);
    for (int i = 1; i < 7; i++){
        val_ana(i-1) = d + 2*a*cos((i*datum::pi)/7);
        for (int j = 1; j < 7; j++){
            //there's gotta be a better way to do this than a nested for loop
            vec_ana(j-1, i-1) = sin((j*i*datum::pi)/7);
        }
    }
    vec_ana = normalise(vec_ana);

    // Print analytical results:
    //std::cout << "Analytical eigenvalues and eigenvectors: "
    //val_ana.print();
    //std::cout << endl;
    //vec_ana.print();
    //std::cout << endl;

    // Print numerical results:
    //std::cout << "Numerical eigenvalues and eigenvectors: "
    //eigval.print();
    //std::cout << endl;
    //eigvec.print();
    //std::cout << endl;

    // Compare numerical method to analytical:
    double tol = 0.0001;
    std::cout << "Comparing numerical (arma::eig_sym) method to analytical:" << endl;
    bool vecsim = check_eigenvectors(eigvec, vec_ana, tol);
    std::cout << "Eigenvectors equivalent True/False: " << vecsim << endl;
    bool valsim = check_eigenvalues(eigval, val_ana, tol);
    std::cout << "Eigenvalues equal True/False: " << valsim << endl;

    // Test max_offdiag_symmetric():
    test_ilode();

    // Now calculate answers using jacobi rotation method:
    double tolerance = 1e-8;
    vec val_jac(N, fill::zeros);
    mat vec_jac(N, N, fill::zeros);
    int maxiter = 40;
    int iterations;
    bool converged;
    jacobi_eigensolver(A, tolerance, val_jac, vec_jac,
                            maxiter, iterations, converged);
    vec_jac = normalise(vec_jac);

    // Print jacobi method results:
    //std::cout << endl;
    //val_jac.print();
    //std::cout << endl;
    //vec_jac.print();

    // Compare jacobi method to analytical:
    std::cout << "Comparing Jacobi method to analytical:" << endl;
    vecsim = check_eigenvectors(vec_jac, vec_ana, tol);
    std::cout << "Eigenvectors equivalent True/False: " << vecsim << endl;
    valsim = check_eigenvalues(val_jac, val_ana, tol);
    std::cout << "Eigenvalues equal True/False: " << valsim << endl;

}
