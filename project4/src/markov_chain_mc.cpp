#include "markov_chain_mc.hpp"

// Compute expectation values
void expectation_values(int L, double T, bool Ordered, int n_cycles, int n_steps,
                        double& eps, double& mag, double& heat, double& sus)
{
  double E2 = 0; // To accumulate total energy**2 exp. val.
  double E = 0; // To accumulate total energy exp.val.
  double M2 = 0; // To accumulate total magn.**2 exp.val.
  double M = 0; // To accumulate total |magn.| exp.val.

  int N = L*L; // Number of spins in the lattice

  // Do 'n_cycles' monte carlo cycles
  for (int j=0; j<n_cycles; j++){
    // Initialize new lattice
    Lattice config(L, T, Ordered);
    // Do one monte carlo cycle, length: 'n_steps'
    for (int i=0; i<n_steps; i++){
      config.spin_flip();
    }
    // Accumulate energy and magnetization expectation values:
    double energy = config.get_energy();
    double magnetization = abs(config.get_magnetization());
    E += energy/n_cycles;
    M += magnetization/n_cycles;
    E2 += energy*energy/n_cycles;
    M2 += magnetization*magnetization/n_cycles;
  }
  // Calculate the values:
  eps = E/N; // Units [J]
  mag = M/N; // Units [1] or [spin]
  heat = (E2 - E*E)/(T*T*N); // units [kB]
  sus = (M2 - M*M)/(T*N); // units [1/J] or [spin^2 / J]
}
