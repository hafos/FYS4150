#include "timecheck.hpp"

void timecheck()
{
  bool check = 0;
  int L = 100;
  int N = L*L;
  int n_cycles = 10;
  double T = 1;
  double energy;
  Lattice config(L, T, check);
  std::clock_t start = std::clock();
  //
  for (int j=0; j<n_cycles; j++){
    // Do one monte carlo cycle, length: 'n_steps'
    for (int i=0; i<N; i++){
      config.spin_flip();
    }
    // Accumulate energy and magnetization expectation values:
    double energy = config.get_energy();
    double magnetization = abs(config.get_magnetization());
    std::cout << energy << " " << magnetization << endl;
    energy = config.total_E_;
    magnetization = abs(config.total_M_);
    std::cout << energy << " " << magnetization << endl << endl;
  }
  //
  std::clock_t end = std::clock();
  double timeused = 1.*(end-start)/CLOCKS_PER_SEC;
  std::cout << "timeused = " << timeused << " seconds " << endl;

  // bool check = 1;
  // int L = 1000;
  // double T = 1;
  // double energy;
  // Lattice config(L, T, check);
  // std::clock_t start = std::clock();
  // for (int i=0; i<10000; i++){
  //   config.spin_flip();
  // }
  // std::clock_t end = std::clock();
  // double timeused = 1.*(end-start)/CLOCKS_PER_SEC;
  // std::cout << "timeused = " << timeused << " seconds " << endl;

  //start = std::clock();
  //for (int i=0; i<100; i++){
  //  energy = config.get_energy_test();
  //}
  //end = std::clock();
  //timeused = 1.*(end-start)/CLOCKS_PER_SEC;
  //std::cout << "timeused = " << timeused << " seconds " << endl;
}
