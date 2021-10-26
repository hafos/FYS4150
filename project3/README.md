Compile with
g++ main.cpp src/Particle.cpp src/PenningTrap.cpp -I include/ -larmadillo -std=c++11 -o main.exe

Manually make the folder 'test_results' before running, or it won't save the files
test_results will contain test result data

The PenningTrap class and Particle class header and source files are
located in src/ and include/

plot_results.py plots frequency data
plot_tests.py plots the test result data
