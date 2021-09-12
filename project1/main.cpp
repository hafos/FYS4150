#include "main.hpp"


int n = 100; // grid points in square matrix
vector<double> a(n-1, -1); // subdiagonal
vector<double> b(n, 2); // diagonal
vector<double> c(n-1, -1); // superdiagonal
vector<double> g(n, 0);

vector<double> x_array(n, 0);
double h = 1.0/(n-1);
double h_squared = pow(h, 2.0); // to lessen FLOPs
double x = 0;


int main(){

    // output file
    std::string filename = "problem7_output.dat";
    std::ofstream ofile;
    ofile.open(filename);
    // spacing in outputfile
    int width = 12; int decimals = 4;

    // filling array
    for (int i = 0; i < n; i++){
    	x_array[i] = x;
        g[i] = f(x_array[i])*h_squared;
    	x += h;
    }

    vector<double> v = general_algorithm(n, a, b, c, g);
    //for (int i = 0; i < n; i++){
    //    cout << v[i] << endl;
    //}

    //  Writing solutions to file
  	for (int i=0; i < n; i++){
  		ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << x_array[i]
  			  << std::setw(width) << std::setprecision(decimals) << std::scientific << v[i]
  			  << std::endl;
  	}

    ofile.close();
    return 0;
}
