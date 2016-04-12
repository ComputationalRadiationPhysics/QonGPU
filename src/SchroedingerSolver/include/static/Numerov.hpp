#ifndef NUMEROV_H_
#define NUMEROV_H_

#include <assert.h>
#include <cmath>
#include <list>
#include <iterator>
#include "hdf5.h"
#include "hdf5_hl.h"

#define DEBUG(x) std::cout<<x<<std::endl;
#define STATUS(x) std::cout<<x<<"...";
#define ENDSTATUS std::cout<<"DONE!"<<std::endl;

#define CHUNKSIZE 100
using namespace std;
__host__ __device__  double V(double x, double t,double z) {
  return 2*z/sqrt(1+x*x);
};
//NumerovKernel to iterate from right to left!
__global__ void iter2(double* psi,
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
    while(tid < nx * ne) {;
        E = dE * (double) tid / (double) nx;
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
//iteration from left to right
__global__ void iter1(double* psi,
					  size_t nx,
					  size_t ne,
					  double xmax,
					  double xmin,
					  double z,
					  double Es,
					  double dE) {

	//get location inside the thread set
	int tid = nx * ( threadIdx.x + blockIdx.x * blockDim.x);
	//define how many steps to jump after each calculation
	int offset = nx * blockDim.x * gridDim.x;
	//define the constants of the simulation
	double dx = fabs(xmax-xmin)/((double)nx);
    double E = Es;
	double heff = 1.0;
	//fn stands for factor n and helps us to structure the code
	double f1, f2, f3;
	//beginning of iterations loop
	while( tid < ne * nx) {

		E += tid * dE / (double)(nx);
		//This loop implements the numerov method
        printf("Initial conditions Psi0 = %f \n",psi[tid]);
		for(auto i = 2; i < nx ;i++){
			f1 = 1.0 / (1.0 + dx * dx  / 12.0 *  (heff * ( V( ( i + 1) * dx + xmin, 0 , z) - E)));
			f2 = ( 1.0 - 5.0 * dx * dx / 12.0 * ( heff *( V( i * dx + xmin, 0, z) - E)));
			f3 = ( 1.0 + dx * dx / 12.0 * ( heff *( V( ( i - 1) * dx + xmin, 0, z) - E)));
			psi[ i + 1 + tid] = f1 * ( 2.0 * psi[ i + tid] * f2 - psi[ i - 1 + tid] * f3);
		};
		tid+=offset;
	};
};

class Numerov: protected StaticSolver1D {

public:
    void solve();
    Numerov(Params1D *pa, std::complex<double> *ps);
    Numerov();
    ~Numerov();

protected:
    vector<double> chunk;
    vector<double> cache;
    list<vector<double>> results;
    list<double> eval;
    void savelevels();
    bool sign(double s);
    void bisect(int j);
    const int nx,ne;
    int z;
    double E;
    const double xmax,xmin;
    Params1D *param;
    
};

#endif 
