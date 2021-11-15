#include "Lattice.hpp"

//constructor
Lattice::Lattice(int L_in)
{
  L_ = L_in;
  //set seed for random number generator using system clock:
  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator_.seed(seed);

  //generate matrix of spins
  std::uniform_int_distribution<int> initial(0,1); //uniform integer distribution

  spin_config_ = imat(L_, L_, fill::none);
  spin_config_.imbue( [&]() { return initial(generator_); } ); //fill spin_config
  spin_config_.replace(0, -1);
}

//return the configuration
imat Lattice::get_config()
{
  return spin_config_;
}

//flip one random spin
void Lattice::spin_flip()
{
  std::uniform_int_distribution<int> row(0, L_-1);
  std::uniform_int_distribution<int> col(0, L_-1);
  int i = row(generator_);
  int j = col(generator_);
  spin_config_(i, j) = -1*spin_config_(i, j);
}
