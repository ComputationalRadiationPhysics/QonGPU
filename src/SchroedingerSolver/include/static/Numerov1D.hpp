#ifndef NUMEROV1D_H
#define NUMEROV1D_H


#include <vector>
#include "math.h"
#include "hdf5_hl.h"

#define DEBUG(x) std::cout<<"Debug "<<x<<endl;


__device__ inline  double V(double x,double t){
	return rsqrt(1+x*x);
};

//This function executes the Numerov method!
__global__ void Numerov(double* psi, int ne,int nx,double  xmax, double xmin){
    
	int tid = nx*(threadIdx.x+blockIdx.x*blockDim.x);
	int offset = nx*blockDim.x*gridDim.x;
	double dx=xmax/((double)nx);
	double E;
	double f1,f2,f3;
	while(tid<nx*ne){
		E=((double)tid)/((double)nx*ne)*(-1.0);
		for(int i = 2; i < nx;i++){
			f1=1/((1+dx*dx/6*(E-V(i*dx,0.0))));
			f2=(1-5/6*dx*dx*(E-V((i-1)*dx,0)));
			f3=(1+dx*dx*(E-V((i-2)*dx,0)));
			psi[i]=f1*(psi[i-1]*f2+psi[i-2]*f3);	
		};
		tid+=offset;
	};
};





class Numerov1D: protected StaticSolver1D {

public:
	void solve();
	void allocresult();
	void bisect();
	Numerov1D(Params1D* pa,std::complex<double>* ps):param(pa),p(ps){
    
	};
	
	Numerov1D(){};
protected:
	std::complex<double>* psil;
	std::complex<double>* p;
	Params1D* param;
	Potential1D pot;
	void tempprint(double* temp,Params1D* p);
};

void Numerov1D::solve(){
	// use local variables for simulation parameters
    int nx = param->getnx();
    int ne = param->getne();
    double xmax = param->getxmax();
    double xmin = param->getxmin();
	// use a new std::vector<double> to save the calculated energy levels in
    std::vector<double> energy();
	//allocate temporal arrays for double precision calculations
    double* temppsi0 = (double*) malloc(sizeof(double)*ne*nx);
    double* dev_psi;
	/*
	 *Use very small initial values since without this initilization
	 *the whole solution would be constant zero(this is because the 
	 *numerov algorithm is a "marching" algorithm) . 
	 *We're doing this in the first steps in order to show what the
	 *for loop only does!
	 */ 
    temppsi0[0]=0;
    temppsi0[1]=1e-10;
	for(int i = nx; i < ne*nx; i+=ne){
		temppsi0[i-1]=0;
		temppsi0[i]=1e-10;
	};
	cudaMalloc((void**)&dev_psi, sizeof(double)*ne*nx);
	cudaMemcpy(dev_psi,temppsi0,sizeof(double)*ne*nx,cudaMemcpyHostToDevice);
	dim3 grid(256);
    dim3 block(3);
	Numerov<<<grid,block>>>(dev_psi,ne,nx,xmax,xmin);
    cudaMemcpy(temppsi0,dev_psi,sizeof(double)*ne*nx,cudaMemcpyDeviceToHost);
	tempprint(temppsi0,param);
    cudaFree((void*) dev_psi);
	free(temppsi0);
};

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
	
	dim = (p->getnx())*(p->getne());
	//print the simulation results in the HDF5 file
	H5LTmake_dataset(fileid, "/numres", 1, &dim, H5T_NATIVE_DOUBLE, tempdata);
	//close the file
	H5Fclose(fileid);
	//free memory
	free(tempdata);
};	
#endif
