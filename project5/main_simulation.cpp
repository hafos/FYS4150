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
  bool save_prob; // Is True if we want to save the probability only (1 file instead of 2)
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
          mysstream >> h_in >> dt_in >> Tmax_in >> xc_in >> sx_in >> px_in >> yc_in >> sy_in >> py_in >> v0_in >> n_slits >> save_prob >> filename;
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
  else // Parameters from command-line
  {
    // Check number of command-line arguments

    if (argc != 14)  // Expect 13 command-line arguments
    {
      // Get the name of the executable file
      std::string executable_name = argv[0];

      std::cerr << "Error: Wrong number of input arguments." << std::endl;
      std::cerr << "Usage: " << executable_name << " <double h> <double dt> <double T> <double xc> <double sx> <double px> <double yc> <double sy> <double py> <double v0> <int n_slits> <bool save_prob> <string filename>" << std::endl;
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
    save_prob = atoi(argv[12]); // If 1 saves only probability. If 0 saves real+imag.
    filename = argv[13]; // filename to save data to (base)
  }

  // Prepare filenames :
  std::string filename_real;
  std::string filename_imag;
  if (save_prob == 0)
  {
    filename_real = filename + "_real.bin";
    filename_imag = filename + "_imag.bin";
  }

  filename += ".bin";


  // Print out the simulation parameters that were read in:
  // This should help prevent mistakes
  std::cout << endl;
  std::cout << "Configuration :" << endl;
  std::cout << "----------------------------------" << endl;
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
  std::cout << "---------------------------------" << endl;
  std::cout << "saving probability only (on/off): " << save_prob << endl;
  if (save_prob)
  {
    std::cout << "prob only filename : " << filename << endl;
  }
  else
  {
    std::cout << "filename_real : " << filename_real << endl;
    std::cout << "filename_imag : " << filename_imag << endl << endl;
  }

  int M = 1/h_in + 1; // Number of elelemts in the box
  int N = M-2; // Number of elements without the borders
  int N2 = N*N;

  // These are a few tests that can help catch unfortunate mistakes in the input:
  if ((M-1)*h_in !=1){
    // Make sure that the step size fits with the simulation box.
    std::cout << "Please specify a step length h so that 1/h is an integer" << endl;
    return 1;
  }

  if (dt_in <= 0){
    // Safeguard from infinitely running loop later
    std::cout << "Please specify a positive, nonzero time step size" << endl;
    return 1;
  }

  if (n_slits < 0){
    // If negative the program would throw out-of-bounds errors later
    std::cout << "Please specify a positive number of slits" << endl;
    return 1;
  }


  Schrodinger syst(M, v0_in, n_slits); // sets up potential and system

  syst.U_init(xc_in, yc_in, sx_in, sy_in, px_in, py_in); // sets up initial state matrix
  cx_mat U = syst.wave_function(); // matrix form
  cx_vec u_vec = mat2vec(M, U); // vector form

  sp_mat V = syst.potential();
  sp_cx_mat A(N2, N2);
  sp_cx_mat B(N2, N2);
  initialize_matrices(M, h_in, dt_in, V, A, B); // sets up Crank-Nicolson

  // 1. store initial condition as first slice of cube
  // 2. Loop over time and add slices to cube
  // 3. Save cube as binary file after ended loop

  double t = 0;
  if (save_prob)
  {
    cube P_n(M, M, 1, fill::zeros);
    P_n.slice(0) = real(U%conj(U));
    while (t < Tmax_in)
    {
      u_vec = compute_next_step(u_vec, A, B);
      U = vec2mat(M, u_vec);
      mat p = real(U%conj(U));
      P_n = join_slices(P_n, p);
      t = t+dt_in;
    }
    // Save as binary data file :
    P_n.save(filename);
  }
  else
  {
    cube Ureal(M, M, 1, fill::zeros);
    cube Uimag(M, M, 1, fill::zeros);
    Ureal.slice(0) = real(U);
    Uimag.slice(0) = imag(U);
    while (t < Tmax_in)
    {
      u_vec = compute_next_step(u_vec, A, B);
      U = vec2mat(M, u_vec);
      Ureal = join_slices(Ureal, real(U));
      Uimag = join_slices(Uimag, imag(U));
      t = t+dt_in;
    }
    // Save as binary data file :
    Ureal.save(filename_real);
    Uimag.save(filename_imag);
  }

}
