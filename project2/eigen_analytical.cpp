//
//  eigen_analytical.cpp
//  
//
//  Created by Katherine Lee on 27/09/2021.
//

#include <stdio.h>
#include <armadillo>

vec eigenvalue(double N, double a, double d){
    vec value(N);
    for (int i = 1; i < N+1; i++){
        value(i-1) = d + 2*a*cos((i*datum::pi)/(N+1));
    }
    return value;
}

mat eigenvector(double N, double a, double d){
    mat vector(N, N);
    for (int i = 1; i < N+1; i++){
        for (int j = 1; j < N+1; j++){
            vector(j-1, i-1) = sin((j*i*datum::pi)/(N+1));
        }
    }
    return vector;
}
