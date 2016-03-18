#ifndef NUMEROV1D_H
#define NUMEROV1D_H



/*
__global__ void Numerov(cuDoubleComplex* psi, cuDoubleComplex ne1,cuDoubleComplex nx1,cuDoubleComplex xm1){
	int nx = (int)nx1.x;
	int xm = (int)xm1;
	int 
	int tid = nx*(threadIdx.x+blockIdx.x*blockDim.x);
	int offset = nx*blockDim.x*gridDim.x;
	cuDoubleComplex dx=xm/((cuDoubleComplex)nx);
	cuDoubleComplex E;
	cuDoubleComplex f1,f2,f3;
	while(tid<nx*ne){
		E=((cuDoubleComplex)tid)/((cuDoubleComplex)nx*ne)*(-1.0);
		for(int i = 2; i < nx;i++){
			f1=1/((1+dx*dx/6*(E-V(i*dx,0.0))));
			f2=(1-5/6*dx*dx*(E-V((i-1)*dx,0)));
			f3=(1+dx*dx*(E-V((i-2)*dx,0)));
			psi[i]=f1*(psi[i-1]*f2+psi[i-2]*f3);	
		};
		tid+=offset;
	};
};
*/
class Numerov: protected StaticSolver1D {

public:
	void solve();
	Numerov(Params1D& p,std::complex<double>* ps){};
	Numerov(){};
};



#endif
