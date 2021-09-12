// Program using the exact soultion of the 1-D Poisson eq. u(x)
//  to solve for x in in range (0, 1) and writing results to a .txt file for plotting

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cmath>

// Exact solution to 1-D poisson
double u(double x){
	return 1 - (1 - exp(-10))*x - exp(-10*x);
}

int main()
{
	// output file
	std::string filename = "problem2_output.dat";
	std::ofstream ofile;
	ofile.open(filename);
	// spacing in outputfile
	int width = 12; int decimals = 4;

	// variables for x array
	int N = 1000;  // number of data points to evaluate
	double xmin = 0; double xmax = 1;
	double h = (xmax - xmin)/N;  // step length

	std::vector<double> x_array(N, 0);  // x-array
	double x = xmin;
	double y = 0;
	// filling array
	for (int i=0; i <= N; i++){
		x_array[i] = x;
		x += h;
		// std::cout << x_array[i] <<"\n"; // to check array
	}
	//  Writing solutions to file
	for (int i=0; i <= N; i++){
		y = u(x_array[i]);

		ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << x_array[i]
			  << std::setw(width) << std::setprecision(decimals) << std::scientific << y
			  << std::endl;
	}
	ofile.close();
	return 0;
}
