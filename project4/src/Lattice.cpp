#include "Lattice.hpp"

//constructor
Lattice::Lattice(int L_in, double T_in)
{
  L_ = L_in;
  arma::vec dE = {8., 4., 0., -4., -8.}; // [J]
  boltzmann_factors_ = arma::exp(-dE/T_in);
  //set seed for random number generator using system clock:
  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator_.seed(seed);

  //generate matrix of spins
  std::uniform_int_distribution<int> initial(0,1); //uniform integer distribution

  spin_config_ = imat(L_+2, L_+2, fill::none);
  spin_config_.imbue( [&]() { return initial(generator_); } ); //fill spin_config
  spin_config_.replace(0, -1);
  // Set borders
  update_borders();
}

//return the configuration
imat Lattice::get_config()
{
  return spin_config_( span(1,L_), span(1,L_) );
}

//flip one random spin
void Lattice::spin_flip()
{
  std::uniform_int_distribution<int> row(1, L_);
  std::uniform_int_distribution<int> col(1, L_);
  int i = row(generator_);
  int j = col(generator_);

  arma::ivec neighbours = {spin_config_(i+1,j), spin_config_(i-1,j),
                                spin_config_(i,j+1), spin_config_(i,j-1)};
  int case_n = arma::sum(arma::abs(spin_config_(i,j) - neighbours)) /2;
  std::uniform_real_distribution<double> U(0.0 ,1.0);
  double r = U(generator_);
  if (r < boltzmann_factors_(case_n)) {
    spin_config_(i, j) = -1*spin_config_(i, j);
    //update_borders(); // Potentially faster than the following if-tests?
    if (i==L_){
      spin_config_(0, j)=spin_config_(i, j);
    }
    else if (i==1){
      spin_config_(L_+1, j)=spin_config_(i, j);
    }
    if (j==L_){
      spin_config_(i, 0)=spin_config_(i, j);
    }
    else if (j==1){
      spin_config_(i, L_+1)=spin_config_(i, j);
    }
  }
}

// Update borders
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
  for (int i=1; i<=L_; i++)
  {
    for (int j=1; j<=L_; j++)
    {
      S += spin_config_(i, j) * spin_config_(i-1, j);
      S += spin_config_(i, j) * spin_config_(i, j-1);
    }
  }
  return -S;
}

// Return total magnetizaton
double Lattice::get_magnetization()
{
  return arma::sum(arma::sum(spin_config_(span(1,L_), span(1,L_))));
}
