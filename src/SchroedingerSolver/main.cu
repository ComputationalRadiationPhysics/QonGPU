// Copyright Maxmilian Böhme

#include <iostream>

#include <complex>

#include "include/SimDef.hpp"

int main(int argc, char** argv) {
  std::complex<double> xma = 50.0;
  std::complex<double> xmi = -50.0;
  std::complex<double> tmi = 0.0;
  std::complex<double> tma = 5.0;
  Params1D p(xma, xmi, tma, tmi, 1e7, 2, 10,1);
  SimDef<Numerov1D, TimeOperator1D, IOHandle, Core1D, Params1D, 1> s(&p);
  s.staticsolve();
  return EXIT_SUCCESS;
}








