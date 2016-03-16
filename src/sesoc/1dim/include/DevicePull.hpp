#ifndef DEVICEPULL
#define DEVICEPULL


/*This class provides only
 *functionality to fetch data
 *from the Device in the CUDA
 *programming model
 *
 */ 


template<class T>
class DevicePull{

public:
	void pull(T& dest,T& src ,size_t count);

};

template<>
class DevicePull<double>{

	void pull(double* dest, double* src,size_t count){

		cudaMemcpy(dest,src,count,cudaMemcpyDeviceToHost);
	};

};
#endif