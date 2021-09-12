#ifdef main_hpp
#define main_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

double f(double x){
    return 100*exp(-10*x);
}

void general_algorithm(vector<double> a, vector<double> b, vector<double> c, vector<double> g, vector<double> v){
    // forward subst.
    double fraction; //  to lessen FLOPs
    g[0] =
    for (int i = 0; i < n; i++){
        fraction = a[i]/b[i];
        b[i+1] = b[i+1] - fraction*c[i];
        g[i+1] = g[i+1] - fraction*g[i];
    }
    // backward subst.
    v[-1] = g[-1]/b[-1];
    // remembering indexing starting with 0 so n-2 starts from second last element
    for (int i = n-2; 0 <= i; i--){
        v[i] = (g[i] - c[i]*v[i+1])/b[i];
    }
    return;
}


#endif
