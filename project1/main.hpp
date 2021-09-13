#ifndef main_hpp
#define main_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;

double f(double x){
    return 100*exp(-10*x);
}

vector<double> general_algorithm(int n, vector<double> a, vector<double> b, vector<double> c, vector<double> g){
    vector<double> v(n, 0);
    // forward subst.
    double fraction; //  to lessen FLOPs
    for (int i = 1; i < n; i++){
        fraction = a[i-1]/b[i-1];
        b[i] = b[i] - c[i-1]*fraction;
        g[i] = g[i] - g[i-1]*fraction;
    }
    // backward subst.
    v[n-1] = g[n-1]/b[n-1];
    // remembering indexing starting with 0 so n-2 starts from second last element
    for (int i = n-2; 0 <= i; i--){
        v[i] = (g[i] - c[i]*v[i+1])/b[i];
    }
    return v;
}

vector<double> specialized_algorithm(int n, vector<double> b, vector<double> g){
    vector<double> v(n, 0);
    // forward subst.
    for (int i = 1; i < n; i++){
        g[i] = g[i] + g[i-1]/b[i-1];
    }
    // backward subst.
    v[n-1] = g[n-1]/b[n-1];
    // remembering indexing starting with 0 so n-2 starts from second last element
    for (int i = n-2; 0 <= i; i--){
        v[i] = (g[i] + v[i+1])/b[i];
    }
    return v;
}


#endif
