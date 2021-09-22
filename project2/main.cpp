//
//  On Mac, compile in the terminal using command:
//  g++ main.cpp -std=c++11 -o main.exe -larmadillo
//

#include <iostream>
#include <armadillo>
#include <cmath>

using namespace arma;


int main() {
    //Problem 3: initialize a tridiagonal 6x6 matrix A
    //with diagonal 2/h^2 and sub/superdiagonal -1/h^2
    //h = 1/(N+1) so here it's 1/7
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
    
    //initialize arrays to contain the eigenvectors and eigenvalues
    vec eigval;
    mat eigvec;
    //solve the eigenvalue problem
    eig_sym(eigval, eigvec, A);
    eigval = normalise(eigval);
    eigvec = normalise(eigvec);
    
    //now check with the analytical result:
    vec val_ana(6);
    mat vec_ana(6, 6);
    for (int i = 1; i < 7; i++){
        val_ana(i-1) = d + 2*a*cos((i*datum::pi)/7);
        for (int j = 1; j < 7; j++){
            //there's gotta be a better way to do this than a nested for loop
            vec_ana(j-1, i-1) = sin((j*i*datum::pi)/7);
        }
    }
    val_ana = normalise(val_ana);
    vec_ana = normalise(vec_ana);
    
    eigvec.print();
    std::cout << endl;
    vec_ana.print(); //whyyyyy are the signs different
    
    std::cout << "eigenvalues equal: " << approx_equal(eigval, val_ana, "absdiff", 0.0001) << endl;
    std::cout << "eigenvectors equal: " << approx_equal(eigvec, vec_ana, "absdiff", 0.0001) << endl;
}
