// #include "main.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

int n = 10; // grid points in square matrix
vector<double> a(n-1, -1); // subdiagonal
vector<double> b(n, 2); // diagonal
vector<double> c(n-1, -1); // superdiagonal
vector<double> g(n, 0);

vector<double> x_array(n, 0);
double h = 1.0/n;
double h_squared = pow(h, 2.0); // to lessen FLOPs
double x = 0;

// double f(double x);
// void general_algorithm(vector<double> a, vector<double> b, vector<double> c, vector<double> g, vector<double> v);

double f(double x){
    return 100*exp(-10*x);
}

vector<double> general_algorithm(vector<double> a, vector<double> b, vector<double> c, vector<double> g){
    vector<double> v(n, 0);
    // forward subst.
    double fraction; //  to lessen FLOPs
    for (int i = 0; i < n; i++){
        fraction = a[i]/b[i];
        b[i+1] = b[i+1] - fraction*c[i];
        g[i+1] = g[i+1] - fraction*g[i];
    }
    // backward subst.
    v[n-1] = g[n-1]/b[n-1];
    // remembering indexing starting with 0 so n-2 starts from second last element
    for (int i = n-2; 0 <= i; i--){
        v[i] = (g[i] - c[i]*v[i+1])/b[i];
    }
    return v;
}

int main(){

    // output file
    // std::string filename = "problem2_output.dat";
    // std::ofstream ofile;
    // ofile.open(filename);
    // spacing in outputfile
    int width = 12; int decimals = 4;

    // filling array
    for (int i = 0; i <= n; i++){
    	x_array[i] = x;
        g[i] = f(x_array[i])/h_squared;
    	x += h;
    }

    vector<double> v = general_algorithm(a, b, c, g);
    for (int i = 0; i < n; i++){
        cout << v[i] << endl;
    }
    // ofile.close();
    return 0;
}
