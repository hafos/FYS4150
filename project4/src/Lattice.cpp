#include "Lattice.hpp"

// Constructor
Lattice::Lattice(int L_in, double T_in, bool Ordered)
{
  L_ = L_in;
  // The 5 energy differences:
  delta_E_ = {8., 4., 0., -4., -8.}; // [J]
  // The 5 corresponding boltzmann factors :
  boltzmann_factors_ = arma::exp(-delta_E_/T_in); // = p(new state) / p(old)
  if (Ordered == 0)
  {
    // Fill spin config with random spins
    // Set seed for random number generator using system clock:
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator_.seed(seed);

    // Generate matrix of spins
    std::uniform_int_distribution<int> initial(0,1); //uniform integer distribution

    spin_config_ = imat(L_+2, L_+2, fill::none);
    spin_config_.imbue( [&]() { return initial(generator_); } ); //fill spin_config
    spin_config_.replace(0, -1);

    // Set outer borders
    update_borders();
  } else
  {
    // Fill spin config with only ones
    spin_config_ = imat(L_+2, L_+2, fill::ones);
  }
  total_E_ = get_energy();
  total_M_ = get_magnetization();
  //std::cout<< total_E_<<endl;
  //std::cout<< total_M_<<endl;
}

// Return the configuration
imat Lattice::get_config()
{
  return spin_config_( span(1,L_), span(1,L_) );
}

// Pick one random spin, flip it if accepted
void Lattice::spin_flip()
{
  // Choose random spin in lattice:
  std::uniform_int_distribution<int> row(1, L_);
  std::uniform_int_distribution<int> col(1, L_);
  int i = row(generator_);
  int j = col(generator_);

  // Find its neighbours and the corresponding case (which dE it is)
  arma::ivec neighbours = {spin_config_(i+1,j), spin_config_(i-1,j),
                                spin_config_(i,j+1), spin_config_(i,j-1)};
  int case_n = arma::sum(arma::abs(spin_config_(i,j) - neighbours)) /2;

  // Generate random r between 0 and 1
  std::uniform_real_distribution<double> U(0.0 ,1.0);
  double r = U(generator_);

  // Accept/reject step:
  if (r <= boltzmann_factors_(case_n)) {
    spin_config_(i, j) = -1*spin_config_(i, j);
    total_E_ += delta_E_(case_n);
    total_M_ += 2*spin_config_(i, j);
    update_borders(); // Potentially faster (with compiler optimization)
    // than the following if-tests?
    // update_borders() seems to be around as fast with -02 and -O3 but slower without
    //if (i==L_){
    //  spin_config_(0, j)=spin_config_(i, j);
    //}
    //else if (i==1){
    //  spin_config_(L_+1, j)=spin_config_(i, j);
    //}
    //if (j==L_){
    //  spin_config_(i, 0)=spin_config_(i, j);
    //}
    //else if (j==1){
    //  spin_config_(i, L_+1)=spin_config_(i, j);
    //}
  }
}

// Update all borders
void Lattice::update_borders()
{
  spin_config_.col(0) = spin_config_.col(L_);
  spin_config_.col(L_+1) = spin_config_.col(1);
  spin_config_.row(0) = spin_config_.row(L_);
  spin_config_.row(L_+1) = spin_config_.row(1);
}

// Return total energy
double Lattice::get_energy()
{
  double S = 0;
  for (int j=1; j<=L_; j++) // Column-major order in armadillo
  {
    arma::ivec tmp_col = spin_config_.col(j); // To reduce memory traffic
    arma::ivec tmp_col_left = spin_config_.col(j-1); // To reduce memory traffic
    for (int i=1; i<=L_; i++) // Do not loop over the outer borders
    {
      S += tmp_col(i) * tmp_col(i-1); // add upward spin interaction
      S += tmp_col(i) * tmp_col_left(i); // and leftward spin interaction
    }
  }
  return -S;
}

// Return total magnetizaton
double Lattice::get_magnetization()
{
  return arma::sum(arma::sum(spin_config_(span(1,L_), span(1,L_))));
}
