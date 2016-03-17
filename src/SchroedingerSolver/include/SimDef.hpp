#ifndef SIMDEF_HPP
#define SIMDEF_HPP
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "boost/static_assert.hpp"
#include "boost/type_traits.hpp"
#include "AllHeader.hpp"
template <class StaticSolver,class TimeOperator,class IOHandler,class Potential,class Params,int dim>
class SimDef
{


	
public:
	SimDef(Params p){};
	void staticSolve();
	void timerev();
	void printres();


	
private:

	StaticSolver s;
	TimeOperator t;
	IOHandler io;
	Potential p;
	Params da;



//Watch for the right cases:
	inline void checktypes(){
	if(dim==1)
	{
		BOOST_STATIC_ASSERT((is_base_of<Paxrams1D,da>::value));
		BOOST_STATIC_ASSERT((is_base_of<StaticSolver1D,s>::value));
		BOOST_STATIC_ASSERT((is_base_of<TimeOperator1D,t>::value));
		BOOST_STATIC_ASSERT((is_base_of<Potential1D,p>::value));
	};
	if(dim==2){
		BOOST_STATIC_ASSERT((is_base_of<Params2D,da>::value));
		BOOST_STATIC_ASSERT((is_base_of<StaticSolver2D,s>::value));
		BOOST_STATIC_ASSERT((is_base_of<TimeOperator2D,t>::value));
		BOOST_STATIC_ASSERT((is_base_of<Potential2D,p>::value));
	};
	
	if(dim==1){
		BOOST_STATIC_ASSERT((is_base_of<Params3D,da>::value));
		BOOST_STATIC_ASSERT((is_base_of<StaticSolver3D,s>::value));
		BOOST_STATIC_ASSERT((is_base_of<TimeOperator3D,t>::value));
		BOOST_STATIC_ASSERT((is_base_of<Potential3D,p>::value));
	};
	};
};



#endif
