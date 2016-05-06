#include "CrankNicholson1D.hpp"
#include "CNKernels.h"



#define DEBUG2(x) std::cout<<x<<std::endl
CrankNicholson1D::CrankNicholson1D(): nx(0),nt(0), E(0) { }

CrankNicholson1D::CrankNicholson1D(Params1D *_p): param(_p),
                                                  nx( _p->getnx()),
                                                  nt( _p->getnt()),
                                                  E( 0.0),
                                                  chunk_h(_p->getnx()),
                                                  chunkl_d( _p->getnx()),
                                                  chunkr_d( _p->getnx()),
                                                  tmax( _p->gettmax()),
                                                  tmin( _p->gettmin()),
                                                  xmax( _p->getxmax()),
                                                  xmin( _p->getxmin()),
                                                  d(_p->getnx()),
                                                  du(_p->getnx()),
                                                  dl(_p->getnx()),
                                                  inital(_p->getnx())
{}
CrankNicholson1D::~CrankNicholson1D() {}

void CrankNicholson1D::rhs_rt() {

    transform_rhs<<<nx,1>>>(raw_pointer_cast(chunkl_d.data()), raw_pointer_cast(chunkr_d.data()), nx, param->getxmax(), param->getxmin(),tau);
}

void CrankNicholson1D::lhs_rt(double x, double t,

                              cuDoubleComplex* d,
                              cuDoubleComplex* du,
                              cuDoubleComplex* dl) { }

void CrankNicholson1D::cusparse_init() {
    status = cusparseCreate( &handle);
    if(status != CUSPARSE_STATUS_SUCCESS) {
        std::cout<<"Error: Critical cuSparse unable to initialize"<<std::endl;
    }
    else {
        std::cout << "cuSparse successfully initialized!" << std::endl;
    }
}

void CrankNicholson1D::cusparse_destr() {
    // This function just destroys the cuSparse context
    // and throws an error if necessary

    status = cusparseDestroy( handle);
    if(status != CUSPARSE_STATUS_SUCCESS) {
        std::cout << "Error: cuSparse wasn't destroyed correctly!" << std::endl;
    }
    else {
        std::cout << "cuSparse was destroyed correctly!" << std::endl;
    }
}


void CrankNicholson1D::time_solve() {

    // This routine is now slightly longer
    
    const double hbar_m = 1.0;
    const double h = (xmax - xmin) / (double) nx;
    const double tau = (tmin - tmax) / (double) nt;
    const double c =  1 / ( h * h) * hbar_m;
    double t = 0;
    cuDoubleComplex* dev_d = raw_pointer_cast(d.data());
    cuDoubleComplex* dev_du = raw_pointer_cast(du.data());
    cuDoubleComplex* dev_dl = raw_pointer_cast(dl.data());
    cuDoubleComplex* dev_rhs = raw_pointer_cast(chunkr_d.data());
    create_const_diag<<<nx ,1>>>( raw_pointer_cast(dl.data()), raw_pointer_cast(du.data()), c, nx);
    cusparse_init();
    for( auto i = 0; i < nt; ++i) {
        t += tau * (double) i;
        rhs_rt();
        // first perform the RHS Matrix multiplication!
        // Then update the non-constant main-diagonal!
        update_diagl<<<nx,1>>>(raw_pointer_cast(d.data()), tau, h, c, xmin, nx, t);
        // right after that, we can call the cusparse Library
        // to write the Solution to the LHS chunk
        status2 = cusparseZgtsv( handle, nx, 1, dev_dl, dev_d, dev_du, dev_rhs, nx);
        std::cout << status2 << std::endl;
        thrust::copy( chunkl_d.begin(), chunkl_d.end(), chunkr_d.begin());
    }
    cusparse_destr();
    MPI_Finalize();
}

void CrankNicholson1D::setstate(const thrust::host_vector <cuDoubleComplex>& v) {
    // copy initial state into memory!
    thrust::copy(v.begin(), v.end(), inital.begin());
    chunkl_d = inital;
    chunkr_d = inital;
}