Compile for tests with :

g++ main_tests.cpp src/*.cpp -I include/ -larmadillo -std=c++14 -o main_tests.exe

Compile for simulations with:

g++ main_simulation.cpp src/*.cpp -I include/ -larmadillo -std=c++14 -o main_simulation.exe

Run for simulations with:

./main_simulation.exe h dt T xc sx px yc sy py v0 n_slits save_prob filename

where h is the position step size, dt is the time step size, T is the simulation time, (xc, yc) is the center of the initial wave packet, (sx, sy) is the width of the initial  wave packet, (px, py) is the momentum of the initial wave packet, and v0 is the height of the potential barrier. n_slits is the number of slits. The full wavefunction is saved (2 files) if save_prob is 0. The probability is saved (1 file) if save_prob is 1. filename (base name, .bin is appended in code) is where the result is saved.

Or use config file:

./main_simulation.exe configurations/config_file_name.dat
