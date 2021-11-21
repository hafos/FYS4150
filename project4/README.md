Compile for tests using "g++ main_tests.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -o main_tests.exe"


Compile for parallelization tests using
- Linux : "g++ main_time_tests.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -fopenmp -o main_time_tests.exe"
- macOS : (UNTESTED) "g++ main_time_tests.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -Xpreprocessor -fopenmp -o main_time_tests.exe"


Compile for final estimations (take long) using
- Linux : "g++ main_parallelized.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -fopenmp -o main_parallelized.exe"
- macOS : (UNTESTED) "g++ main_parallelized.cpp src/*.cpp -I include/ -larmadillo -O2 -std=c++11 -Xpreprocessor -fopenmp -o main_parallelized.exe"


The Lattice class header and source files are in /include and /src respectively.
