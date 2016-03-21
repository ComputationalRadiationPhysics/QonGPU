#ifndef SIMDEF_HPP
#define SIMDEF_HPP
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "cuComplex.h"
#include "boost/static_assert.hpp"
#include "boost/type_traits/is_base_of.hpp"
#include "thrust/complex.h"
#include "AllHeader.hpp"



template <class StatSolver,class TimeOp,class IO,class Pot,class Para,int dim>
class SimDef{
public:
	SimDef(Params1D *p):da(p),s(p,psi0){
		 
	};
	SimDef(Params2D *p):da(p){
	};
	SimDef(Params3D *p):da(p){
	};
	void staticsolve(){
		std::cout << "Debug 1" << std::endl;
		s.solve();
		DEBUG("CALL 2")
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
	double* npsi;



};


#endif
