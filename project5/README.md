## Repository with code to simulate the 2D time-dependent schrodinger equation

- Compile for basic tests with :

g++ main_tests.cpp src/*.cpp -I include/ -larmadillo -std=c++14 -o main_tests.exe

*Note: main_tests.cpp was mainly a diagnostic for us to see that seperate parts of the code behaved as we expected. Thus it is not referenced in the paper.

- Compile for simulations with:

g++ main_simulation.cpp src/*.cpp -I include/ -larmadillo -std=c++14 -o main_simulation.exe

- Run for simulations with a config file:

./main_simulation.exe configurations/config_file_name.dat

- Run for simulations with parameters in the command line:

./main_simulation.exe h dt T xc sx px yc sy py v0 n_slits save_prob filename

where h is the position step size, dt is the time step size, T is the simulation time, (xc, yc) is the center of the initial wave packet, (sx, sy) is the width of the initial  wave packet, (px, py) is the momentum of the initial wave packet, and v0 is the height of the potential barrier. n_slits is the number of slits. The full wave function is saved (2 files) if save_prob is 0. The probability is saved (1 file) if save_prob is 1. filename (base name, .bin is appended in code) is where the result is saved.



After running all the config_file scenarios spanning problems 7-9, the results can be plotted using plot_simulation.py, which will also produce animations of the 5 scenarios. If one just wants to see the animations they are in the repo in the directory /animations.

To make a file containing two animations that play side by side, use the following command:

ffmpeg -i left_filename.mp4 -i right_filename.mp4 -filter_complex hstack output_filename.mp4
