#include <iostream>
#include "include/DevicePusher.hpp"
#include "include/DeviceAlloc.hpp"
#include "include/DevicePull.hpp"
#include "DataInit.hpp"
#include "Definition.hpp"
#include "Potential.hpp"
#include "NumerovKernel.hpp"
#include "StaticSolver.hpp"
int main(int argc, char **argv){

	SimData* s=Allocation(NE,NX);
	s->energy[0]=10000;

	std::cout << "Everything worked!" <<std::endl;
	std::cout << s->energy[0] << std::endl;
	Destruction(s);
	return 0;
}

