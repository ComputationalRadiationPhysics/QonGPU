#ifndef NUMEROV1D_H
#define NUMEROV1D_H


#include <vector>

#include "math.h"

#include <cmath>

#include "hdf5.h"

#include "hdf5_hl.h"

#include "stdio.h"

#include <assert.h>

#define DEBUG(x) std::cout<<x<<endl;

#define FORWARDS

__device__  double V(double x, double t,double z) {
  return z/sqrt(1+x*x);
};

__global__ void NumerovKernel2(double* psi,
                               int nx,
                               int ne,
                               double xmax,
                               double xmin,
                               double z) {
    //Should always show to the psi_n(x=0) at the energy level E_n
    int tid = nx * ( threadIdx.x + blockIdx.x * blockDim.x );
    int offset = nx * blockDim.x * gridDim.x;
    double dx = fabs( xmax - xmin ) / ( (double) nx );
    double E = -1.0;
    double dE = E / ( (double) ne );
    double f1, f2,f3;
    double heff = 1 / 12 * dx * dx;
    while(tid < nx * ne) {
        E = dE * (double) tid / (double) nx;
        printf("Using the energy: %lf \n",E);
        for (auto i = nx - 1; i >= 2; i -= 1) {
            f1 = 1.0 + heff * (E + V(xmax - (i - 2) * dx , 0, z));
            f2 = 1.0 - 5.0 * heff * (E + V(xmax -( i - 1) * dx, 0, z));
            psi[i - 2 + tid] = (2 * psi[i - 1 + tid] * f2 - psi[i  + tid] * f3) / f1;
        }
        tid += offset;
    }
}

__global__ void NumerovKernel(double* psi,
                              size_t nx,
                              size_t ne,
                              double xmax,
                              double xmin,
                              double z) {
	unsigned int tid = nx*(threadIdx.x+blockIdx.x*blockDim.x);//Should always show to the psi_n(x=0) at the energy level E_n
    uint offset = nx*blockDim.x*gridDim.x;
	double dx = fabs(xmax-xmin)/((double)nx);
	double E = -V(0,0,z);
	double dE=E/((double)ne);
	double f1,f2,f3;
	double heff=1.0;
	while(tid<ne*nx){
		E=tid*dE/(nx);
		for(auto i = 2; i<nx ;i++){
			f1=1/(1+dx*dx*(heff*(E+V((i+1)*dx+xmin,0,z)))/12);
			f2=(1-5*dx*dx/12*(heff*(E+V(i*dx+xmin,0,z))));
			f3=(1+dx*dx/12*(heff*(E+V((i-1)*dx+xmin,0,z))));	
			psi[i+1+tid]=f1*(2*psi[i+tid]*f2-psi[i-1+tid]*f3);
		};
		tid+=offset;
    };
};


class Numerov1D: protected StaticSolver1D {

public:
    void solve();
    void allocresult();
    void bisect();
    Numerov1D(Params1D* pa,std::complex<double>* ps):param(pa),p(ps)
    ,pot(),eindex(100),out(pa){
	nx = (param->getnx());
	ne =  (param->getne());
	psi=(double*) malloc(sizeof(double)*nx*ne);
#ifdef BACKWARDS
	if(psi!=NULL){
	    for(size_t i = nx; i < nx*ne; i+=nx){
		    psi[i]=0;
		    psi[i-1]=-1e-10;
	    };
	}
#endif
#ifdef FORWARDS
    if(psi!=NULL){
	    for(size_t i = 0; i < nx*ne; i+=nx){
	    	psi[i]=0;
		    psi[i+1]=-1e-10;
	    };
	}
#endif
	else{
	    std::cout<<"Error: Allocation failed!"<<std::endl;
	    assert(psi!=NULL);
	}
    };
    Numerov1D(){};
    ~Numerov1D(){
	free(psi);
    };
protected:
    std::complex<double>* psil;
    std::complex<double>* p;
    //initilize the energy index vector
    std::vector<int> eindex;
    double* psi;
    Params1D* param;
    Potential1D pot;
    void tempprint(double* temp,Params1D* p);
    size_t nx,ne;
    bool sign(double x);
    int hindex=0;
    double z=param->getz();
    IOHandle1D out;
};


/*
 *This function provides the Cuda-Kernel call for the Numerov method.
 *It's primarily here to be used as an interface, which helps to allocate memory
 *and call the important functions! 
 */ 
void Numerov1D::solve(){
    double* dev_ptr;
    //Allocate needed Memomry, in this case an nx*ne double array!
    cudaMalloc((void**)&dev_ptr,sizeof(double)*nx*ne);
    //copy Memory to the device
    cudaMemcpy(dev_ptr,psi,sizeof(double)*nx*ne,cudaMemcpyHostToDevice);
    //call the actual Kernel
#ifdef BACKWARDS
    DEBUG("Running backwards iteration")
    NumerovKernel2<<<256,3>>>(dev_ptr,nx,ne,param->getxmax(),param->getxmin(),z);
#endif
#ifdef FORWARDS
    DEBUG("Running forwards iteration")
    NumerovKernel<<<256,3>>>(dev_ptr,nx,ne,param->getxmax(),param->getxmin(),z);
#endif
    //fetch the results into the psi array!
    cudaMemcpy(psi,dev_ptr,sizeof(double)*nx*ne,cudaMemcpyDeviceToHost);
    //free the cuda Memory
    cudaFree(&dev_ptr);
    //use bisection method
#ifdef  FORWARDS
    bisect();
#endif
    //generate output
    tempprint(psi,param);
    
};


/*
 *This function provides temporariy the possibility to print out the
 *wave function! It will late be removed. This has to be here, because the
 *static results should be tested firtst!
 */
void Numerov1D::tempprint(double* temp,Params1D* p){
    hid_t fileid;
    hsize_t dim=7;
    //Create temporal array for parameter location:
    double* tempdata = (double*) malloc(sizeof(double)*7);
    tempdata[0]=p->getxmax();
    tempdata[1]=p->getxmin();
    tempdata[2]=p->gettmax();
    tempdata[3]=p->gettmin();
    tempdata[4]=(double) p->getnx();
    tempdata[5]=(double) p->getnt();
    tempdata[6]=(double) p->getne();
    //create the hdf5 file:
    fileid = H5Fcreate("static_results.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    //Print the parameters in the file
    H5LTmake_dataset(fileid,"/params1d",1,&dim, H5T_NATIVE_DOUBLE,tempdata);
    free(tempdata);
    //rewrite dims
    dim = (p->getnx())*(p->getne());
    //print the simulation results in the HDF5 file
    H5LTmake_dataset(fileid, "/numres", 1, &dim, H5T_NATIVE_DOUBLE,temp);
    //now write the energy indices to the file
    tempdata =(double*) malloc(sizeof(double)*(hindex));
    for(auto i=0;i<hindex;i++){
	tempdata[i]=(double)eindex[i];
	}
	//create the index file
    dim=hindex;
    H5LTmake_dataset(fileid,"/levels",1,&dim,H5T_NATIVE_DOUBLE,tempdata);
    free(tempdata);
    //close the file
    H5Fclose(fileid);
    //free memory
    DEBUG("Finished saving!")	
};

/*
 *Find out if in the searched interval any energy levels were found
 *and include the index as well as the energy in a std::vector<double>.
 *
 */
bool Numerov1D::sign(double x){
    return std::signbit(x);
}
void Numerov1D::bisect(){
    bool act=false;
    bool prev=false;
    double fine=1e-3;
    hindex=0;
    for(auto i=1;i<ne;i++) {
	    act=sign(psi[nx*i+nx-1]);
	    prev=sign(psi[nx*(i-1)+nx-1]);
	    if(!(act==prev)){
	        if(fabs(psi[nx*i+nx-1])< fine) {
                DEBUG((psi[nx * i + nx - 1]))
                eindex[hindex] = i - 1;
                hindex += 1;
	        }
	    else {
		DEBUG("Psi was too large")
	    }
 	}		    		  
    };
  
    DEBUG("----------------------------------")
	DEBUG("Finished!")
	DEBUG(hindex<<" Energy levels were found!")
};

#endif
