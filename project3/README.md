Compile with
g++ main.cpp src/Particle.cpp src/PenningTrap.cpp -I include/ -larmadillo -std=c++11 -o main.exe


RK4_contents_note.txt is my attempt at RK4, was not finished and I prioritized testing euler to
make sure the code structure as a whole was working well. Will fix that later.

euler_test.dat contains test data of euler method. plots.py plots the x and y values.
Two particles are included, not counting interactions. I don't think the results are reasonable yet.
