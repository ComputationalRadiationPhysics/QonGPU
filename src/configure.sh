#!/usr/bin/bash

if [$1 = "-m"];then
    echo"Module loading enabled!"
    module purge
    module load gcc/4.9.2
    module load cuda/7.5
    module load hdf5
    module load cmake/3.3.0
    module load boost
fi


cmake .
mv CMakeFiles/ build/
mv CMakeCache.txt build/
mv cmake_install.cmake build/

echo "Running qsolve"
build/Solver/qsolve



