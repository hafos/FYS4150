#include "markov_chain_mc.hpp"

// Compute expectation values
void expectation_values(int L, double T, bool Ordered, int n_cycles,
                        double& eps, double& mag, double& heat, double& sus)
{
  double E2 = 0; // To accumulate total energy**2 exp. val.
  double E = 0; // To accumulate total energy exp.val.
  double M2 = 0; // To accumulate total magn.**2 exp.val.
  double M = 0; // To accumulate total |magn.| exp.val.

  int N = L*L; // Number of spins in the lattice

  // Initialize random lattice
  Lattice config(L, T, Ordered);

  // Do 'n_cycles' monte carlo cycles
  for (int j=0; j<n_cycles; j++){
    // Do one monte carlo cycle, length: 'n_steps'
    for (int i=0; i<N; i++){
      config.spin_flip();
    }
    // Accumulate energy and magnetization expectation values:
    double energy = config.total_E_; //config.get_energy();
    double magnetization = abs(config.total_M_); //abs(config.get_magnetization());
    E += energy/n_cycles;
    M += magnetization/n_cycles;
    E2 += energy*energy/n_cycles;
    M2 += magnetization*magnetization/n_cycles;

  }
  // Calculate and store the values:
  eps = E/N; // Units [J]
  mag = M/N; // Units [1] or [spin]
  heat = (E2 - E*E)/(T*T*N); // units [kB]
  sus = (M2 - M*M)/(T*N); // units [1/J] or [spin^2 / J]
}

// Compute energy  and magnetization expectation value for each cycle
void energy_cycles(int L, double T, bool Ordered, int n_cycles, vec&e, vec&m)
{
  int N = L*L; // Number of spins in the lattice
  vec energy(n_cycles);
  vec magnetization(n_cycles);
  // Initialize new lattice
  Lattice config(L, T, Ordered);
  // Store intials
  e(0) = config.total_E_/N;
  m(0) = config.total_M_/N;
  // Do 'n_cycles' monte carlo cycles, store energy after each
  for (int j = 0; j<n_cycles; j++){
    // Do one monte carlo cycle, length: 'n_steps'
    for (int i = 0; i<N; i++){
      config.spin_flip();
    }
    //update energy
    double energ = config.get_energy(); // no integer div
    energy(j) = energ/N;
    //update magnetization
    double magnet = abs(config.get_magnetization()); // no integer div
    magnetization(j) = magnet/N;
  }
  //cumulative sum of energies, divided by number of cycles
  energy = cumsum(energy)/linspace(1, n_cycles+1, n_cycles);
  //cumulative sum of magnetizations, divided by number of cycles
  magnetization = cumsum(magnetization)/linspace(1, n_cycles+1, n_cycles);
  for (int i=0; i<energy.n_elem; i++)
  {
    e(i+1) = energy(i);
    m(i+1) = magnetization(i);
  }
}
