#include "timecheck.hpp"

void timecheck()
{
  bool check = 1;
  int L = 1000;
  double T = 1;
  double energy;
  Lattice config(L, T, check);
  std::clock_t start = std::clock();
  for (int i=0; i<10000; i++){
    config.spin_flip();
  }
  std::clock_t end = std::clock();
  double timeused = 1.*(end-start)/CLOCKS_PER_SEC;
  std::cout << "timeused = " << timeused << " seconds " << endl;

  //start = std::clock();
  //for (int i=0; i<100; i++){
  //  energy = config.get_energy_test();
  //}
  //end = std::clock();
  //timeused = 1.*(end-start)/CLOCKS_PER_SEC;
  //std::cout << "timeused = " << timeused << " seconds " << endl;
}
