#ifndef STATICSOLVER
#define STATICSOLVER
#include "DeviceAlloc.hpp"
#include "DevicePusher.hpp"
#include "DevicePull.hpp"
#include "NumerovKernel.hpp"

void SolveForInitial(size_t nx,size_t ne,SimData* data){
	double* dev_psi;
	
	for(int i=0;i<nx*ne;i+=nx){
		data->psi[i]=0;
		data->psi[i+1]=1e-10;
	};
	DeviceAlloc<double> al();
	DevicePusher<double> push();
	DevicePull<double> pull();
	al.alloc(dev_psi,nx*ne);
	
};



#endif
