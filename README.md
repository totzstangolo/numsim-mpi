# numsim-mpi
mkdir -p build 
cd build
rm -rf CMakeFiles/ CMakeCache.txt cmake_install.cmake Makefile
cmake DEBUG_VISU=ON ..
make
time mpirun -n 4 ./NumSim

rm -rf CMakeFiles/ CMakeCache.txt cmake_install.cmake Makefile && cmake DEBUG_VISU=ON .. && make && time mpirun -n 4 ./numsim
