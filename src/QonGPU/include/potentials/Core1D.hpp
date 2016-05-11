#ifndef CORE1D_H_
#define CORE1D_H_

#include "Potential1D.hpp"

class Core1D:Potential1D {
public:
    // This is just a simple functor to give us the Value of the 
    // Core Potential!
    __host__ __device__ inline  double operator()(double x, double t) {
	return 1/sqrt(x*x+1);
    }

};
#endif
