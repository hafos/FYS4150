//test particle count
double sim_time = 1; // [microseconds]
double B0_in = 9.65*1e1;
double V0_in = 9.65*1e8;
double d_in = 1e4;

bool Interactions = 1;
PenningTrap PT(B0_in, V0_in, d_in, Interactions);

double q_in = 1;
double m_in = 1;
vec r_in = { 1, 0, 0 };
vec v_in = { 0, 0.1, 0 };
vec r_in2 = {0, 1, 0};
vec v_in2 = { 0, 0, 0 };

Particle p_in = Particle(q_in, m_in, r_in, v_in);
PT.add_particle(p_in);
Particle p_in2 = Particle(q_in, m_in, r_in2, v_in2);
PT.add_particle(p_in2);

double t = 0;
int count = 0;

int maxiter = 10;
while (t <= sim_time)
{
  PT.evolve_RK4(dt);
  std::cout << PT.count_particles() << endl;
  t += dt;
  count +=1;
  if (count >= maxiter){
    std::cout << "Maxiter, t:" << t << endl;
    break;
  }
}

test random particles:
double q_in = 1;
double m_in = 1;
double B0_in = 9.65*1e1;
double V0_in = 9.65*1e8;
double d_in = 1e4;
int n = 10;
bool Interactions = 1;

PenningTrap PT(B0_in, V0_in, d_in, Interactions);
PT.add_n_particles(n, q_in, m_in);
std::cout << PT.count_particles() << endl;
