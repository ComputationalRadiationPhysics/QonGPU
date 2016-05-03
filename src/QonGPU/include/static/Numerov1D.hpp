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
  return 2*z/sqrt(1+x*x);
};
//NumerovKernel to iterate from right to left!
__global__ void NumerovKernel2(double* psi,
                               int nx,
                               int ne,
                               double xmax,
                               double xmin,
                               double z) {
    //Should always show to the psi_n(x=0) at the energy level E_n
    int tid = nx * ( threadIdx.x + blockIdx.x * blockDim.x );
    int offset = nx * blockDim.x * gridDim.x;
    double dx = ( xmax - xmin ) / nx;
    double E = -V(0,0,z);
    double dE = E/ne;
    double f1, f2,f3;
    double heff = dx * dx / 12.0;
    while(tid < nx * ne) {
        E = dE * tid /  nx;
        for(auto i = nx-1; i>1;i--) {
            //refine the first pre-factor for psi_(n-1)
            f1 = 1.0 + heff * (E + V(xmax - (i-2) * dx, 0, z));
            //redefine the next prefactor for psi_n
            f2 = 1.0 - 5.0 * heff * (E + V(xmax - (i-1) *dx, 0, z));
            //redefine at lat the third prefactor for psi_(n+1)
            f3 = 1.0 + heff * (E + V(xmax - i * dx , 0, z));
            //finally let's calculate the new value psi_n-1=1/f1 * (2 * psi_n * f2 - psi_n+1 * f3)
            psi[ i + tid - 2 ] = 1 / f1 * (2 * psi[i+tid-1] * f2 - psi[ i+ tid ] * f3);
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
	int tid = nx*(threadIdx.x+blockIdx.x*blockDim.x);
    //Should always show to the psi_n(x=0) at the energy level E_n
	int offset = nx*blockDim.x*gridDim.x;
	double dx = fabs(xmax-xmin)/((double)nx);
	double E = V(0,0,z);
	double dE=E/((double)ne);
	double f1,f2,f3;
	double heff=1.0;
	while(tid<ne*nx){
		E=tid*dE/(nx);
		for(auto i = 2; i<nx ;i++){
			f1=1.0/(1.0+dx*dx*(heff*(V((i+1)*dx+xmin,0,z)-E))/12.0);
			f2=(1.0-5.0*dx*dx/12.0*(heff*(V(i*dx+xmin,0,z)-E)));
			f3=(1.0+dx*dx/12.0*(heff*(V((i-1)*dx+xmin,0,z)-E)));
			psi[i+1+tid]=f1*(2.0*psi[i+tid]*f2-psi[i-1+tid]*f3);
		};
		tid+=offset;
    };
};


class Numerov1D: protected StaticSolver1D {

public:
    void solve();
    void bisect();
    void bisect2();
    Numerov1D(Params1D* pa) : param(pa)
    , pot(), eindex(100), out(pa)
    {
	nx = (param->getnx());
	ne =  (param->getne());
	psi=(double*) malloc(sizeof(double)*nx*ne);
#ifdef BACKWARDS
	if(psi!=NULL){
	    for(size_t i = nx-1; i < nx*ne; i+=nx){
		    psi[i]=0;
		    psi[i-1]=1e-10;
	    };
	}
#endif
#ifdef FORWARDS
    if(psi!=NULL){
	    for(size_t i = 0; i < nx*ne; i+=nx){
	    	psi[i]=0;
		    psi[i+1]=1e-10;
	    };
	}
#endif
    else{
	        std::cout<<"Error: Allocation failed!"<<std::endl;
            std::cout<<"Please make sure to have at least: "<<nx*ne*4/1024/1024<<"megabytes RAM availible!"<<endl;
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
    size_t nx,ne;
    bool sign(double x);
    int hindex=0;
    double z=param->getz();
    IOHandle1D out;
    void tempprint(double* temp,Params1D* p);

};


/*
 *This function provides the Cuda-Kernel call for the Numerov method.
 *It's primarily here to be used as an interface, which helps to allocate memory
 *and call the important functions! 
 */ 

void Numerov1D::solve(){
    double* dev_ptr;
    //Allocate required device Memomry, in this case an nx*ne double array!
    cudaMalloc((void**)&dev_ptr,sizeof(double)*nx*ne);
    //copy Memory to the device
    cudaMemcpy(dev_ptr,psi,sizeof(double)*nx*ne,cudaMemcpyHostToDevice);
    //call the actual Kernel
#ifdef BACKWARDS
    DEBUG("Running backwards iteration!")
    NumerovKernel2<<<256,4>>>(dev_ptr,nx,ne,param->getxmax(),param->getxmin(),z);
#endif
#ifdef FORWARDS
    DEBUG("Running forwards iteration!")
    NumerovKernel<<<256,4>>>(dev_ptr,nx,ne,param->getxmax(),param->getxmin(),z);
#endif
    //fetch the results into the psi array!
    cudaMemcpy(psi,dev_ptr,sizeof(double)*nx*ne,cudaMemcpyDeviceToHost);
    //free the cuda Memory
    cudaFree(&dev_ptr);
    //use bisection method
#ifdef  FORWARDS
    bisect();
#endif
#ifdef BACKWARDS
    bisect2();
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

void Numerov1D::bisect2() {
    auto act = false;
    auto prev = false;
    double fine = 1e-4;
    hindex=0;
    for(auto i = 1;i < ne; i++){
        act = sign(psi[nx*i]);
        prev = sign(psi[nx*(i-1)]);
        if(!(act==prev)) {
            if (fabs(psi[nx * i]) <= fine) {
                DEBUG(psi[nx * i])
                eindex[hindex] = i - 1;
                hindex += 1;
            }

            else {
                DEBUG("Psi was too large")
            }
        }
    }
    DEBUG("----------------------------------")
    DEBUG("Finished!")
    DEBUG(hindex<<" Energy levels were found!")
}

#endif
