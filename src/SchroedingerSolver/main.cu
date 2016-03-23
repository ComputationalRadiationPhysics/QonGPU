
#include <iostream>
#include "include/SimDef.hpp"
#include <complex>


int main(int argc,char** argv){
	
	std::complex<double> xma=20.0;
	std::complex<double> xmi=0.0;
	std::complex<double> tmi=0.0;
	std::complex<double> tma=5.0;
	Params1D p(xma,xmi,tma,tmi,1e6,2,50);
	DEBUG("CALL 1")
	SimDef<Numerov1D,TimeOperator1D,IOHandle,Potential1D,Params1D,1> s(&p);
	s.staticsolve();
	DEBUG("CALL 3")
	return 0;
}








