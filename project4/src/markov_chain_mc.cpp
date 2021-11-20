#include "markov_chain_mc.hpp"

// Compute expectation values
void expectation_values(int L, double T, bool Ordered, int n_cycles, int n_burnin,
                        double& eps, double& mag, double& heat, double& sus)
{
  double E2 = 0; // To accumulate total energy**2 exp. val.
  double E = 0; // To accumulate total energy exp.val.
  double M2 = 0; // To accumulate total magn.**2 exp.val.
  double M = 0; // To accumulate total |magn.| exp.val.

  int N = L*L; // Number of spins in the lattice

  // Initialize random lattice
  Lattice config(L, T, Ordered);

  // Burn-in
  for (int j = 0; j<n_burnin; j++){
    for (int k = 0; k<N; k++){
      config.spin_flip();
    }
  }

  // Do 'n_cycles' monte carlo cycles
  for (int j=0; j<n_cycles; j++){
    // Do one monte carlo cycle, length N
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

// Get energy and magnetization samples
void individual_samples(int L, double T, bool Ordered, int n_cycles, int n_burnin,
                        vec&e, vec&m)
{
  int N = L*L; // Number of spins in the lattice
  // Initialize new lattice
  Lattice config(L, T, Ordered);
  // Store intials
  e(0) = config.total_E_/N; // Not really interesting if we remove burn-in
  m(0) = config.total_M_/N;
  // Burn-in
  for (int j = 0; j<n_burnin; j++){
    for (int k = 0; k<N; k++){
      config.spin_flip();
    }
  }
  // Do 'n_cycles' monte carlo cycles, store energy after each
  for (int j = 0; j<n_cycles; j++){
    // Do one monte carlo cycle, length N
    for (int i = 0; i<N; i++){
      config.spin_flip();
    }
    // Update energy
    double energ = config.get_energy(); // no integer div
    e(j+1) = energ/N;
    // Update magnetization
    double magnet = abs(config.get_magnetization()); // no integer div
    m(j+1) = magnet/N;
  }
}
