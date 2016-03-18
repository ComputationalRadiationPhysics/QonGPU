#include <iostream>
#include "SimDef.hpp"
#include <complex>

int main(int argc,char** argv){
	
	std::complex<double> xma=5.0;
	std::complex<double> xmi=0.0;
	std::complex<double> tmi=0.0;
	std::complex<double> tma=5.0;
    Params1D p(xma,xmi,tma,tmi,0,0,0);
	SimDef<Numerov,TimeOperator1D,IOHandle,Potential1D,Params1D,1> s(&p);
	
	return 0;
}








