//
// Created by max on 22/04/16.
//

#include "CrankNicholson1D.hpp"
#include "CNKernels.h"

CrankNicholson1D::CrankNicholson1D(): nx(0),nt(0), E(0) { }

CrankNicholson1D::CrankNicholson1D(Params1D *_p, vector <cuDoubleComplex> _v): param( _p),
                                                                               nx( _p->getnx()),
                                                                               nt( _p->getne()),
                                                                               E( 0.0),
                                                                               chunkl_d( _p->getnx()),
                                                                               chunkr_d( _p->getnx()),
                                                                               chunk_h( _v),
                                                                               tmax( _p->gettmax()),
                                                                               tmin( _p->gettmin()),
                                                                               xmax( _p->getxmax()),
                                                                               xmin( _p->getxmin()),
                                                                               d(nx),
                                                                               du(nx),
                                                                               dl(nx)
{
    thrust::copy( chunk_h.begin(), chunk_h.end(), chunkr_d.begin());
}
CrankNicholson1D::~CrankNicholson1D() { }

void CrankNicholson1D::rhs_rt() {

    transform_rhs<<<nx,1>>>(raw_pointer_cast(chunkl_d.data()), raw_pointer_cast(chunkr_d.data()), nx, param->getxmax(), param->getxmin(),tau);
}

void CrankNicholson1D::lhs_rt(double x, double t,
                              cuDoubleComplex* d,
                              cuDoubleComplex* du,
                              cuDoubleComplex* dl) {

}

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
    const double c =  tau / 2 / ( h * h) * hbar_m;


    double t = 0;
    cuDoubleComplex* dev_d = raw_pointer_cast(d.data());
    cuDoubleComplex* dev_du = raw_pointer_cast(du.data());
    cuDoubleComplex* dev_dl = raw_pointer_cast(dl.data());
    cuDoubleComplex* dev_rhs = raw_pointer_cast(chunkr_d.data());

    create_const_diag<<<nx ,1>>>( raw_pointer_cast(dl.data()), raw_pointer_cast(du.data()), c, nx);
    cusparse_init();
    for( auto i = 0; i < nt; ++i) {
        t += tau * ( double) i;
        // first perform the RHS Matrix multiplication!
        rhs_rt();
        // Then update the non-constant main-diagonal!
        update_diagl<<<nx,1>>>(raw_pointer_cast(d.data()), tau, h, c, xmin, nx, t);
        // right after that, we can call the cusparse Library
        // to write the Solution to the LHS chunk
        status2 = cusparseZgtsv( handle, nx, 1, dev_dl, dev_d, dev_du, dev_rhs, nx);
        if(status2 != CUSPARSE_STATUS_SUCCESS) {
            std::cout << status2 << std::endl;
        }
        thrust::copy( chunkl_d.begin(), chunkl_d.end(), chunkr_d.begin())Q;
    }
    cusparse_destr();

}