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

void CrankNicholson1D::rhs_rt( double x, double tau) {

    transform_rhs<<<nx,1>>>(raw_pointer_cast(chunkl_d.data()), raw_pointer_cast(chunkr_d.data()), nx, param->getxmax(), param->getxmin(),tau);
}

void CrankNicholson1D::lhs_rt(double x, double t,
                              cuDoubleComplex* d,
                              cuDoubleComplex* du,
                              cuDoubleComplex* dl) {

}

void CrankNicholson1D::time_solve() {
    // This routine is now slightly longer

    const double hbar_m = 1.0;
    const double h = (xmax - xmin) / (double) nx;
    const double tau = (tmin - tmax) / (double) nt;
    const double c =  tau / 2 / ( h * h) * hbar_m;

    create_const_diag<<<nx ,1>>>( raw_pointer_cast(dl.data()), raw_pointer_cast(du.data()), c, nx);


}