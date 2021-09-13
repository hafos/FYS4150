#include "main.hpp"

vector<vector<double>> calculate_the_values(int n, bool timing=false){
  vector<double> a(n-1, -1); // subdiagonal
  vector<double> b(n, 2); // diagonal
  vector<double> c(n-1, -1); // superdiagonal
  vector<double> g(n, 0);

  vector<double> x_array(n, 0);
  double h = 1.0/(n+1);
  double h_squared = pow(h, 2.0); // to lessen FLOPs
  double x = h;

  // filling array
  for (int i = 0; i < n; i++){
    x_array[i] = x;
      g[i] = f(x_array[i])*h_squared;
    x += h;
  }

  if (timing == true){
      if (n == 10){
          // vector<vector<double>> time;
          std::string filename2 = "timing_run.dat";
          std::ofstream otimefile;
          otimefile.open(filename2);
          otimefile << "n: " << n << endl;
          for (int i = 0; i < 10; i++){
              clock_t t1 = clock();
              vector<double> gen_time = general_algorithm(n, a, b, c, g);
              clock_t t2 = clock();
              double gen_time_elapsed = ((double) (t2 - t1))/CLOCKS_PER_SEC;
              // time[0].push_back(gen_time_elapsed);

              clock_t t3 = clock();
              vector<double> special_time = specialized_algorithm(n, b, g);
              clock_t t4 = clock();
              double special_time_elapsed = ((double) (t4 - t3))/CLOCKS_PER_SEC;
              // time[0].push_back(special_time_elapsed);

              otimefile << gen_time_elapsed << " "
              << special_time_elapsed
              << std::endl;
          }
          otimefile.close();
      }
      else{
           // vector<vector<double>> time;
           std::ofstream otimefile;
           otimefile.open("timing_run.dat", std::ofstream::app);
           otimefile << "n: " << n << endl;
           for (int i = 0; i < 10; i++){
               clock_t t1 = clock();
               vector<double> gen_time = general_algorithm(n, a, b, c, g);
               clock_t t2 = clock();
               double gen_time_elapsed = ((double) (t2 - t1))/CLOCKS_PER_SEC;
               // time[0].push_back(gen_time_elapsed);

               clock_t t3 = clock();
               vector<double> special_time = specialized_algorithm(n, b, g);
               clock_t t4 = clock();
               double special_time_elapsed = ((double) (t4 - t3))/CLOCKS_PER_SEC;
               // time[0].push_back(special_time_elapsed);

               otimefile << gen_time_elapsed << " "
               << special_time_elapsed
               << std::endl;
           }
           otimefile.close();
       }
      }
  vector<double> v = general_algorithm(n, a, b, c, g);
  vector<vector<double>> xv(2,vector<double> (n, 0));
  for (int i = 0; i < n; i++){
      xv[0][i] = x_array[i];
      xv[1][i] = v[i];
      //    cout << v[i] << endl;
  }
  return xv;
}



int main(){
    // output file
    std::string filename = "problem7_output.dat";
    std::ofstream ofile;
    ofile.open(filename);
    // spacing in outputfile
    int width = 21; int decimals = 9;

    //  Writing solutions to file
    int n = 10; // grid points in square matrix
    vector<vector<double>> xv = calculate_the_values(n, true);
  	for (int i=0; i < n; i++){
  		ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[0][i]
  			  << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[1][i]
  			  << std::endl;
    }
    ofile << std::endl;

    n = 100; // grid points in square matrix
    xv = calculate_the_values(n, true);
  	for (int i=0; i < n; i++){
  		ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[0][i]
  			  << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[1][i]
  			  << std::endl;
  	}
    ofile << std::endl;

    n = 1000; // grid points in square matrix
    xv = calculate_the_values(n, true);
  	for (int i=0; i < n; i++){
  		ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[0][i]
  			  << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[1][i]
  			  << std::endl;
    }
    ofile << std::endl;

    n = 10000; // grid points in square matrix
    xv = calculate_the_values(n, true);
  	for (int i=0; i < n; i++){
  		ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[0][i]
  			  << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[1][i]
  			  << std::endl;
    }
    ofile << std::endl;
    n = 100000; // grid points in square matrix
    xv = calculate_the_values(n, true);
  	for (int i=0; i < n; i++){
  		ofile << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[0][i]
  			  << std::setw(width) << std::setprecision(decimals) << std::scientific << xv[1][i]
  			  << std::endl;
    }
    //ofile << std::endl;

    ofile.close();
    return 0;
}
