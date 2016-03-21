#ifndef NUMEROV1D_H
#define NUMEROV1D_H


#include <vector>
#include "math.h"
#include "hdf5_hl.h"
#include "stdio.h"
#define DEBUG(x) std::cout<<"Debug "<<x<<endl;


__device__ double V(double x,double t){
	return 0.5*x*x;
};


/*
//This function executes the Numerov method!
__global__ void Numerov(double* psi, int ne,int nx,double  xmax, double xmin){
    
	int tid = nx*(threadIdx.x+blockIdx.x*blockDim.x);
	int offset = nx*blockDim.x*gridDim.x;
	double dx=xmax/((double)nx);
	double E=0;
	double f1,f2,f3;
	while(tid<nx*ne){
		E=((double)tid)/((double)nx*ne)*(-1.0);
		for(int i = 2; i < nx;i++){
			f1=1/((1+dx*dx/6*(E-V(((double)i)*dx,0.0))));
			f2=(1-5/6*dx*dx*(E-V(((double)i-1)*dx,0)));
			f3=(1+dx*dx*(E-V(((double)i-2)*dx,0)));
			psi[i+tid]=f1*(psi[i-1+tid]*f2+psi[i-2+tid]*f3);	
		};
		tid+=offset;
	};
};
*/

__global__ void NumerovKernel(double* psi, int nx,int ne, double xmax,double xmin){
	int tid = nx*(threadIdx.x+blockIdx.x*blockDim.x);
	int offset = nx*blockDim.x*gridDim.x;
	double dx = xmax/((double)nx);
	double E;
	double f1,f2,f3;
	double xpos;
	while(tid<ne*nx){
		E=1/((double)ne)*((double)tid/nx);
		for(int i = 2; i<nx-1;i++){
			xpos=dx*(double)i;
			f1=1/(1+dx*dx/6*(E-V(xpos,0)));
			f2=(1-5.0/6.0*dx*dx*(E-V(xpos+dx,0)));
			f3=(1+dx*dx/6.0*(E-V(xpos-dx,0)));
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
	Numerov1D(Params1D* pa,std::complex<double>* ps):param(pa),p(ps){
		 nx = (param->getnx());
		 ne =  (param->getne());
		psi=(double*) malloc(sizeof(double)*nx*ne);
		for(int i = 0; i < nx*ne; i+=nx){
			psi[i]=0;
			psi[i+1]=1e-5;
		};
	};
	Numerov1D(){};
	~Numerov1D(){
		free(psi);
	};
protected:
	std::complex<double>* psil;
	std::complex<double>* p;
	std::vector<double> energy;
	double* psi;
	Params1D* param;
	Potential1D pot;
	void tempprint(double* temp,Params1D* p);
	int nx,ne;
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
	NumerovKernel<<<256,1>>>(dev_ptr,nx,ne,param->getxmax(),param->getxmin());
	//fetch the results into the psi array!
	cudaMemcpy(psi,dev_ptr,sizeof(double)*nx*ne,cudaMemcpyDeviceToHost);
	//free the cuda Memory
	cudaFree(&dev_ptr);
	//generate output
	tempprint(psi,param);
	//use bisection method
	bisect();
};


/*
 *This function provides temporariy the possibility to print out the
 *wave function! It will late be removed. This has to be here, because the
 *results should be tested firtst!
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
	//rewrite dims
	dim = (p->getnx())*(p->getne());
	//print the simulation results in the HDF5 file
	H5LTmake_dataset(fileid, "/numres", 1, &dim, H5T_NATIVE_DOUBLE,temp);
	//close the file
	H5Fclose(fileid);
	//free memory
	DEBUG("OUTPUT_CALL")
	free(tempdata);
};

/*
 *Find out if in the searched interval any energy levels were found
 *and include the index as well as the energy in a std::vector<double>
 *
 */
/*
void Numerov1D::bisect(){
	int pre=signbit(psi[nx]);
	int preindex=nx;
	for(int i = 2*nx; i<ne*nx;i+=ne){
		if( signbit(psi[i])!= pre){
			pre=signbit(psi[i]);
			preindex=i;
			std::cout<<"No Level on this index pair!"<<std::endl;
			std::cout<<"Sign of Psi(pre): "<<pre<<std::endl;
		}
		else{
			if(fabs(psi[preindex])<fabs(psi[i])){
				std::cout<<"Energy Level found!"<<std::endl;
				
			}
			else{
				std::cout<<"Energy Level found!"<<std::endl;
			}
		}
	}	
	};
*/
#endif
