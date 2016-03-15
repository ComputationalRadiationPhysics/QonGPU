#ifndef POTENTIAL
#define POTENTIAL

 
__device__ double V(double x,double t){

	return rsqrt(1+x*x); 

};
#endif 
