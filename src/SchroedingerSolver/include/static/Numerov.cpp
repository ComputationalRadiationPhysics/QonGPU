#include "Numerov.hpp"

Numerov::Numerov(Params1D *pa,std::complex<double>* ps): params(pa) {
    nx = params->getnx();
    ne = params->getne();
    z = params->getz();
    psi1=(double*) malloc(sizeof(double)*nx*ne);
    psi2=(double*) malloc(sizeof(double)*nx*ne);
    if(psi1==NULL||psi2==NULL){
        DEBUG("Error! Parameters too big!"<<endl
	      <<"Memory allocation failed!")
        if(psi1==NULL) assert(psi1!=NULL);
	if(psi2==NULL) assert(psi2!=NULL);
    }
    for(auto i = 0; i < nx*ne;i+=nx){
	psi2[i+nx-1]=0;
	psi2[i+nx-2]=1e-10;
	psi1[i]=0;
	psi1[i+1]=1e-10;
    }
}
    
Numerov::Numerov() {}
Numerov::~Numerov() {
    free(psi1);
    free(psi2);
}
void Numerov::solve(){   
    // This function runs a rightward iteration first
    // Therefore we  need to reserver graphics memory
    double* dev_psi1;
    double* dev_psi2;
    cudaMalloc((void**)&dev_psi1,sizeof(double)*nx*ne);
    STATUS("Running inverse iteration...")
    cudaMemcpy(dev_psi1,psi1,sizeof(double)*nx*ne,cudaMemcpyHostToDevice);
    iter1<<<256,4>>>(dev_psi1,nx,ne,0,params->getxmin(),z);
    ENDSTATUS("DONE!")
    cudaMemcpy(psi1,dev_psi1,sizeof(double)*nx*ne,cudaMemcpyDeviceToHost);
    
    cudaFree(dev_psi1);
    cudaMalloc((void**)&dev_psi2,sizeof(double)*nx*ne);
    cudaMemcpy(dev_psi2,psi2,sizeof(double)*nx*ne,cudaMemcpyHostToDevice);
    STATUS("Running regular iteration...")
    iter2<<<256,4>>>(dev_psi2,nx,ne,params->getxmax(),0,z);
    STATUS("DONE!")
    cudaMemcpy(psi2,dev_psi2,sizeof(double)*nx*ne,cudaMemcpyDeviceToHost);
    cudaFree(dev_psi2);
}

bool Numerov::sign(double x){
    return std::signbit(x);
}

void Numerov::bisect() {
    /*
    //show numer of found energy levels
    int hind1=0;
    int hind2=0; 
    //index locations of energy levels
    std::vector<int> en;
    //create boolean expressions for easier statements
    bool exp1;
    bool exp2;
    for(auto i = 1; i < ne; i++) {
        exp1 = (sign(psi2[i*nx])==sign(psi[(i-1)*nx]));
	if(!exp1){
	}
		     
    }*/    
}

void Numerov::savelevels(std::vector<int> en,int hind) {
/*
    //This functions saves all the required parameters
    //and values(only detected energy levels)
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
    // now preparing the required data!
*/
}

