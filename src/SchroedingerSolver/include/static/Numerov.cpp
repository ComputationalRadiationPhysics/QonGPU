#include "Numerov.hpp"
Numerov::Numerov():nx(0),ne(0),xmax(0),xmin(0){}

Numerov::Numerov(Params1D *pa,complex<double>* ps): param(pa),
          chunk(pa->getnx()*CHUNKSIZE),cache(pa->getnx()*CHUNKSIZE),
						    nx(pa->getnx()),
						    ne(pa->getne()),
						    z(pa->getz()),
						    xmax(pa->gettmax()),
						    xmin(pa->getxmin())						    
{
  // usin the init cache as long as necessary
  for(auto it=cache.begin(); it!=cache.end();it+=nx) {
      *(it+nx)=0;
      *(it+nx-1)=1e-10;
  }
  //get the minimal energy E_n is in [-V(0,0,z),0]
  E=-V(0,0,z);
  
}

Numerov::~Numerov(){}

void Numerov::solve(){
    //create a double device ponter
    double* dev_ptr;
    //now we enter a loop
    //which computes a "chunk" of Solutions
    //which gets analyzed and the overwritten
    //after each step
    STATUS("Allocation of graphics memory")
    cudaMalloc((void**)&dev_ptr,sizeof(double)*nx*CHUNKSIZE);
    ENDSTATUS
    double dE = E/(double) ne;
    auto E_lok=0.0;
    for(auto j = 0; j < ne;j+=CHUNKSIZE) {
	//push the chunk vector (with initial conditions
	//into the allocated device memory	
	E_lok+=dE*(double)j;
	cudaMemcpy(dev_ptr,chunk.data(),sizeof(double)*CHUNKSIZE*nx,cudaMemcpyHostToDevice);
	STATUS("Calculating the "<<j<<"-th Chunk")
        iter1<<<256,4>>>(dev_ptr,nx,CHUNKSIZE,xmax,xmin,z,E_lok,dE);
	cudaMemcpy(cache.data(),dev_ptr,sizeof(double)*CHUNKSIZE*nx,cudaMemcpyDeviceToHost);
	ENDSTATUS
	STATUS("Running bisection")
        bisect(j);
        ENDSTATUS
    }
    cudaFree(dev_ptr);
}
void Numerov::savelevels(){}
bool Numerov::sign(double s){
    return std::signbit(s);
}
void Numerov::bisect(int j){
    for(auto it=cache.begin();it!=cache.end();it+=nx){
	if(sign(*(it+nx-1))!=sign(*(it-1))){
	        DEBUG("Signchange detected!")
		if(fabs(*(it+nx-1))<fabs(*(it-1))&&
					 fabs(*(it+nx-1))<1e-3){
		       vector<double> v(it,(it+nx-1));
		       results.push_back(v);
		       eval.push_back((double)j*E/(double)ne);
		       DEBUG("Energy level found at"<<(eval.back()))
		}
	        else if(fabs(*(it-1))<1e-3){
		    vector<double> v(it-nx-1,it-1);
		    results.push_back(v);
		    eval.push_back((double)(j-1)*E/(double)ne);
		    DEBUG("Energy level found at"<<(eval.back()))
		}
	}
	else{
	    DEBUG("No sign change detected")
	    }
    }
			   
}

