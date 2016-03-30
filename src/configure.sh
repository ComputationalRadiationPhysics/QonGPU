#!/usr/bin/bash

cmake .
mv CMakeFiles/ build/
mv CMakeCache.txt build/
mv cmake_install.cmake build/

echo "Running qsolve"
build/Solver/qsolve



