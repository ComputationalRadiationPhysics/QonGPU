#ifndef SIMDEF_HPP
#define SIMDEF_HPP
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cuComplex.h"
#include "boost/static_assert.hpp"
#include "boost/type_traits/is_base_of.hpp"
#include "thrust/complex.h"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "AllHeader.hpp"



template <class StatSolver,class TimeOp,class IO,class Pot,class Para,int dim>
class SimDef
{


	
public:
	SimDef(Params1D *p){
		da=p;
	};
	SimDef(Params2D *p){
	};
	SimDef(Params3D *p){
	};
	void staticsolve(){
    };
	void timerev(){

	};
	void printres(){

	};


	
private:

	StatSolver s;
	TimeOp t;
	IO io;
	Pot p;
	Para *da;
	std::complex<double>* psi0;
	std::complex<double>* psi;
	double npsi;



};


#endif
