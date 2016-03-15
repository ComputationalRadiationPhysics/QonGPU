#ifndef DEVICEALLOC
#define DEVICEALLOC


/*This is a general Wrapper for cudaMalloc
 *included in the Cuda programming model for
 *general types. It's specialized for double 
 *types since these will be mostly used in the
 *Simulation. 
 *
 */ 


template<class T>
class DeviceAlloc{
public:
	void alloc(T& dev_p,size_t size);


};

template<>
class DeviceAlloc<double>{
	void alloc(double* dev_p,size_t size){

		cudaMalloc((void**) dev_p,size);

	};


};
#endif
