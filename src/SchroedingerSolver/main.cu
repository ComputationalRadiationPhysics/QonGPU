// Copyright Maxmilian Böhme

#define DEBUG(x) std::cout<<x<<std::endl;
#define STATUS(x) std::cout<<x<<"...";
#define ENDSTATUS std::cout<<"DONE!"<<std::endl;

#include "include/SimDef.hpp"

int main(int argc, char** argv) {
  std::complex<double> xma = 5.0;
  std::complex<double> xmi = 0.0;
  std::complex<double> tmi = 0.0;
  std::complex<double> tma = 10.0;
  Params1D p(xma, xmi, tma, tmi, 1e3, 2, 1e3,1);
  SimDef<Numerov, TimeOperator1D, IOHandle, Core1D, Params1D, 1> s(&p);
  s.staticsolve();
  return EXIT_SUCCESS;
}








