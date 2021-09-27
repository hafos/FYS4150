//
//  create_tridiag.cpp
//  
//
//  Created by Katherine Lee on 27/09/2021.
//

#include <stdio.h>
#include <armadillo>

mat create_symmetric_tridiagonal(int n, double a, double d){
    mat A = mat(n, n, fill::eye);
    A(0, 0) = d;
    A(0, 1) = a;
    for (int i = 1; i < n-1; i++){
        A(i, i-1) = a;
        A(i, i) = d;
        A(i, i+1) = a;
    }
    A(n-1, n-2) = a;
    A(n-1, n-1) = d;
    
    return A;
}
