//
//  On Mac, compile in the terminal using command:
//  g++ main.cpp -std=c++11 -o main.exe -larmadillo
//

#include <iostream>
#include <armadillo>
#include <cmath>
#include "jacobi.cpp" // includes ilode.cpp
#include "vector_check.cpp"
#include "create_tridiag.cpp"
#include "eigen_analytical.cpp"
#include "scaling.cpp"

using namespace arma;


int main() {
    // Problem 3:
    // initialize a tridiagonal 6x6 matrix A
    // with diagonal 2/h^2 and sub/superdiagonal -1/h^2
    // N = n-1
    // h = 1/n = 1/(N+1)
    int N = 6;
    double h = 1./(N+1);
    double a = -1.0/pow(h, 2);
    double d = 2.0/pow(h, 2);

    mat A = create_symmetric_tridiagonal(N, a, d);

    // Solve the eigenvalue problem numerically:
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);
    eigvec = normalise(eigvec);

    // Compute analytical result:
    vec val_ana = eigenvalue(N, a, d);
    mat vec_ana = eigenvector(N, a, d);
    vec_ana = normalise(vec_ana);

    // Compare numerical method to analytical:
    double tol = 1e-5; // Tolerance of the relative difference
    std::cout << "Problem 3:" << endl;
    std::cout << "Comparing numerical (arma::eig_sym) method to analytical:" << endl;
    std::cout << " (Considered equal if the relative difference is smaller than " << tol << ")" << endl;
    bool vecsim = check_eigenvectors(eigvec, vec_ana, tol);
    std::cout << " Eigenvectors equivalent (True=1/False=0): " << vecsim << endl;
    bool valsim = check_eigenvalues(eigval, val_ana, tol);
    std::cout << " Eigenvalues equal (True=1/False=0): " << valsim << endl;

    // Print analytical results:
    //std::cout << "Analytical eigenvalues and eigenvectors: " << endl;
    //val_ana.print();
    //std::cout << endl;
    //vec_ana.print();
    //std::cout << endl;

    // Print numerical results:
    //std::cout << "Numerical eigenvalues and eigenvectors: " << endl;
    //eigval.print();
    //std::cout << endl;
    //eigvec.print();
    //std::cout << endl;

    // Problem 4
    // Test max_offdiag_symmetric():
    std::cout << endl << "Problem 4:" << endl;
    test_ilode();
    std::cout << endl;

    // Problem 5
    std::cout << "Problem 5: " << endl;
    // Now calculate answers using jacobi rotation method:
    double tolerance = 1e-8; // Absolute tolerance of off-diagonal elements in eigenvalue-matrix
    vec val_jac(N, fill::zeros);
    mat vec_jac(N, N, fill::zeros);
    int maxiter = 40;
    int iterations;
    bool converged;
    jacobi_eigensolver(A, tolerance, val_jac, vec_jac,
                            maxiter, iterations, converged);
    if (converged==1){
      std::cout << "Jacobi rotation method converged after " << iterations;
      std::cout << " iterations. (tolerance = " << tolerance << " )" << endl;
    }
    vec_jac = normalise(vec_jac);

    // Compare jacobi method to analytical:
    std::cout << "Comparing Jacobi method to analytical:" << endl;
    std::cout << " (Considered equal if the relative difference is smaller than " << tol << ")" << endl;
    vecsim = check_eigenvectors(vec_jac, vec_ana, tol);
    std::cout << " Eigenvectors equivalent (True=1/False=0): " << vecsim << endl;
    valsim = check_eigenvalues(val_jac, val_ana, tol);
    std::cout << " Eigenvalues equal (True=1/False=0): " << valsim << endl;

    // Print jacobi method results:
    //std::cout << "Jacobi method results: " << endl;
    //std::cout << endl;
    //val_jac.print();
    //std::cout << endl;
    //vec_jac.print();

    // Problem 6 a)
    /* commented out so not to generate anew for every run
       uncomment to run when you want the .dat file for generating plot
       with problem6_plot.py
    */

    // std::string filename = "problem6_a_output.dat";
    // std::ofstream ofile;
    // ofile.open(filename);
    // int iteration;
    // for (int i=2; i < 121; i+=2){
    //     scaling_N(i, iteration);
    //     ofile << "N : " << i << std::endl;
    //     ofile << "iterations : "
    //           << iteration << std::endl;
    // }
    // ofile.close();

    // Problem 6 b)

    // std::string filename = "problem6_b_output.dat";
    // std::ofstream ofile;
    // ofile.open(filename);
    // int iteration;
    // for (int i=2; i < 123; i+=4){
    //     ofile << "N : " << i << std::endl;
    //     for (int j=0; j < 11; j++){
    //         scaling_dense(i, iteration);
    //         ofile << "iterations : "
    //         << iteration << std::endl;
    //     }
    // }
    // ofile.close();

    // Problem 7
    // first solve for n = 10
    double h_10 = 1./10.;
    double a_10 = -1./pow(h_10, 2);
    double d_10 = 2./pow(h_10, 2);
    mat A_10 = create_symmetric_tridiagonal(9, a_10, d_10);

    // solve eigenvalue problem
    vec val_10(9, fill::zeros);
    mat vec_10(9, 9, fill::zeros);
    jacobi_eigensolver(A_10, tolerance, val_10, vec_10,
                       125, iterations, converged);

    //normalise and sort
    val_10 = normalise(val_10);
    uvec indices = sort_index(val_10);
    val_10 = sort(val_10);

    vec_10 = normalise(vec_10);
    vec_10 = vec_10.cols(indices(span::all));

    // Save results:
    vec_10.save("vec_10.bin");
    val_10.save("val_10.bin");

    // compute analytically for comparison:
    vec val_ana_10 = eigenvalue(9, a_10, d_10);
    val_ana_10 = normalise(val_ana_10);
    uvec indices_ana = sort_index(val_ana_10);
    val_ana_10 = sort(val_ana_10);

    mat vec_ana_10 = eigenvector(9, a_10, d_10);
    vec_ana_10 = normalise(vec_ana_10);
    vec_ana_10 = vec_ana_10.cols(indices_ana(span::all));

    val_ana_10.save("val_ana_10.bin");
    vec_ana_10.save("vec_ana_10.bin");


    // Repeat for n = 100
    double h_100 = 1./100.;
    double a_100 = -1./pow(h_100, 2);
    double d_100 = 2./pow(h_100, 2);
    mat A_100 = create_symmetric_tridiagonal(99, a_100, d_100);

    vec val_100(99, fill::zeros);
    mat vec_100(99, 99, fill::zeros);
    jacobi_eigensolver(A_100, tolerance, val_100, vec_100,
                       15000, iterations, converged);

    //normalise and sort
    val_100 = normalise(val_100);
    uvec indices_100 = sort_index(val_100);
    val_100 = sort(val_100);

    vec_100 = normalise(vec_100);
    vec_100 = vec_100.cols(indices_100(span::all));

    // Save results:
    vec_100.save("vec_100.bin");
    val_100.save("val_100.bin");

    // compute analytically for comparison:
    vec val_ana_100 = eigenvalue(99, a_100, d_100);
    val_ana_100 = normalise(val_ana_100);
    mat vec_ana_100 = eigenvector(99, a_100, d_100);
    vec_ana_100 = normalise(vec_ana_100);
    val_ana_100.save("val_ana_100.bin");
    vec_ana_100.save("vec_ana_100.bin");
}
