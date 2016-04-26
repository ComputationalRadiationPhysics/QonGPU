//
// Created by max on 22/04/16.
//

#include "CrankNicholson1D.hpp"

__device__ __host__ cuDoubleComplex inline pot(double x){
    return make_cuDoubleComplex(1/sqrt(x*x+1),0);
}

__global__ void transfrm_rhs(cuDoubleComplex* rhs,size_t nx,double x) {

    int ind = threadIdx.x + blockDim.x * blockId.x;
    int oset = blockDim.x * gridDim.x;
    cuDoubleComplex x1, x2, x3;
    while(ind < nx) {
        
    }
}

CrankNicholson1D::CrankNicholson1D(): nx(0),nt(0), E(0) { }

CrankNicholson1D::CrankNicholson1D(Params1D *_p, vector <cuDoubleComplex> _v): param(_p),
                                                                               nx(_p->getnx()),
                                                                               nt(_p->getne()),
                                                                               E(0.0),
                                                                               chunk_d(_p->getnx()),
                                                                               chunk_h(_p->getne()){}
CrankNicholson1D::~CrankNicholson1D() { }

void CrankNicholson1D::rhs_rt() {

}