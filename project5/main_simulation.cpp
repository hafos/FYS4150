#include <iostream> // Printing to terminal and file
#include <fstream> // Open and read file
#include <sstream> // stringstream
#include <string> // strings

#include "CrankNicolson.hpp"
#include "Schrodinger.hpp"

using namespace arma;

int main(int argc, char* argv[])
{
  // Get configuration from config file or command-line:
  double h_in, dt_in, Tmax_in, xc_in, sx_in, px_in, yc_in, sy_in, py_in, v0_in;
  int n_slits;
  std::string filename;
  if (argc == 2)
  {
    // Read in from config file
    std::string config_file = argv[1];

    std::fstream myfile;
    myfile.open(config_file);
    if (myfile.is_open()) // Make sure that the file was opened correctly
    {
      std::string line;
      int linecount = 0; // To make sure there is only one line with numbers ..

      while (std::getline(myfile, line))
      {
        if (line.at(0) == '#') // Skip lines that start with #
        {
          continue;
        }
        else
        {
          if (linecount !=0)
          {
            std::cout << "Warning: More than one configuration in config file " << std::endl;
          }
          // Save the contents of the line :
          std::stringstream mysstream(line);
          mysstream >> h_in >> dt_in >> Tmax_in >> xc_in >> sx_in >> px_in >> yc_in >> sy_in >> py_in >> v0_in >> n_slits >> filename;
          linecount += 1;
        }
      }
    }
    else
    {
      std::cerr << "Unable to open the file : " << config_file << std::endl;
      return 1;
    }
    myfile.close();
  }
  else
  {
    // Check number of command-line arguments

    if (argc != 13)  // Expect 12 command-line arguments
    {
      // Get the name of the executable file
      std::string executable_name = argv[0];

      std::cerr << "Error: Wrong number of input arguments." << std::endl;
      std::cerr << "Usage: " << executable_name << " <double h> <double dt> <double T> <double xc> <double sx> <double px> <double yc> <double sy> <double py> <double v0> <int n_slits> <string filename>" << std::endl;
      std::cerr << "Or     " << executable_name << " <string config_filename>" << std::endl;

      // Exit program with non-zero return code to indicate a problem
      return 1;
    }

    // read in parameters from command line
    h_in = atof(argv[1]); // position step size (not Planck constant)
    dt_in = atof(argv[2]); // time step size
    Tmax_in = atof(argv[3]); // end of time interval
    xc_in = atof(argv[4]); // center of initial wave packet (x)
    sx_in = atof(argv[5]); // initial width of wave packet (x)
    px_in = atof(argv[6]); // initial momentum of wave packet (x)
    yc_in = atof(argv[7]); // center of initial wave packet (y)
    sy_in = atof(argv[8]); // initial width of wave packet (y)
    py_in = atof(argv[9]); // initial momentum of wave packet (y)
    v0_in = atof(argv[10]); // potential barrier height
    n_slits = atoi(argv[11]); // number of slits, must be positive or 0
    filename = argv[12]; // filename to save data to
  }

  std::cout << endl;
  std::cout << "Configuration :" << endl;
  std::cout << "---------------------------------------" << endl;
  std::cout << "       h : " << h_in << endl;
  std::cout << "      dt : " << dt_in << endl;
  std::cout << "       T : " << Tmax_in << endl;
  std::cout << "      xc : " << xc_in << endl;
  std::cout << "      sx : " << sx_in << endl;
  std::cout << "      px : " << px_in << endl;
  std::cout << "      yc : " << yc_in << endl;
  std::cout << "      sy : " << sy_in << endl;
  std::cout << "      py : " << py_in << endl;
  std::cout << "      v0 : " << v0_in << endl;
  std::cout << " n_slits : " << n_slits << endl;
  std::cout << "filename : " << filename << endl << endl;

  int M = 1/h_in + 1;
  int N = M-2;
  int N2 = N*N;
  Schrodinger syst(M, v0_in, n_slits); // sets up potential and system

  syst.U_init(xc_in, yc_in, sx_in, sy_in, px_in, py_in); // sets up initial state matrix
  cx_mat U = syst.wave_function(); // matrix form
  //std::cout << U << endl;
  cx_vec u_vec = mat2vec(M, U); // vector form

  sp_mat V = syst.potential();
  sp_cx_mat A(N2, N2);
  sp_cx_mat B(N2, N2);
  initialize_matrices(M, h_in, dt_in, V, A, B); // sets up Crank-Nicolson

  // loop over time
  // store initial condition as first slice of probability density cube
  cube P_n(M, M, 1, fill::zeros);
  P_n.slice(0) = real(U%conj(U));

  double t = 0;
  mat p(M, M);
  while (t < Tmax_in)
  {
    u_vec = compute_next_step(u_vec, A, B);
    U = vec2mat(M, u_vec);
    p = real(U%conj(U));
    P_n = join_slices(P_n, p);
    t = t+dt_in;
  }

  //std::cout << P_n.tail_slices(1) << endl;

  //P_n.each_slice( [](mat& X){ std::cout << accu(X) << endl; } );

  // save as binary data file
  P_n.save(filename);


}