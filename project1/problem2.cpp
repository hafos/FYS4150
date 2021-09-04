// Plot this exact solution of one-dimensional 
// poisson equation

#include <iostream>
#include <vector>

int main()
{
	int N = 10;
	std::vector<double> x(N,0);
	for (int i=0; i<N; i++){
		x[i] = i;
		std::cout << x[i] <<" ";
	}
}
