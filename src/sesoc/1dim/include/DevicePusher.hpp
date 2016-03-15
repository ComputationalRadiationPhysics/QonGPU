#ifndef DEVICEPUSHER
#define DEVICEPUSHER

template<class T>
class DevicePusher{

public:
	DevicePusher();
    void copy(T& dest,T& src,size_t count);
private:
	~DevicePusher();
};


/*
 *This Template specialization is only temporary I will
 *specify it later on. 
 *
 */ 

template<>
class DevicePusher<double>{

	void copy(double* dest,double* src,size_t count){
		cudaMemcpy(dest,src,count,cudaMemcpyHostToDevice);

	};


};


#endif
