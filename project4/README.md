Compile 'main_tests.exe' for tests using: "g++ main_tests.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -o main_tests.exe"


Compile 'main_time_tests.exe' for parallelization tests using:
- Linux : "g++ main_time_tests.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -fopenmp -o main_time_tests.exe"
- macOS : "g++ main_time_tests.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -Xpreprocessor -fopenmp -o main_time_tests.exe -lomp"


Compile 'main_parallelized.exe' for final estimations using:
- Linux : "g++ main_parallelized.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -fopenmp -o main_parallelized.exe"
- macOS : "g++ main_parallelized.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -Xpreprocessor -fopenmp -o main_parallelized.exe -lomp"


The Lattice class header and source files are in /include and /src respectively. It contains methods to initialize a lattice, flip spins according to the criteria, and return it's properties.

The markov_chain_mc .cpp and .hpp files in /include and /src contain the methods to sample states and compute the expectation values, using the Lattice class.

To plot the test data from main_tests.exe, we use plot_tests_w_analytical.py.

To plot the result data from main_parallelized.exe we use plot_results.py.
