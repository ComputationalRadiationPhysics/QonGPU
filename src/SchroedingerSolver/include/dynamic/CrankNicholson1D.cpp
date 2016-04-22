//
// Created by max on 22/04/16.
//

#include "CrankNicholson1D.hpp"

__host__ __device__ double Pot(double x,double t){
    return 2 / sqrt(1 + x * x);
}
__global__ void PrepRHS(double* v) {
    // This function preperates the Right hand side of the CN algorithm

}

__global__ void PrepLHS(double* d1, double* d2, double* d3, size_t nx, double t) {
    // This function preperates the Left hand side of the CN algorithm
}

CrankNicholson1D::CrankNicholson1D() : nx(0),
                                       nt(0) {

}

CrankNicholson1D::CrankNicholson1D(Params1D& p):
                                    nx(p.getnx()),
                                    nt(p.getnt()),
                                    params(&p),
                                    frames(nx*TCHUNK){

}

CrankNicholson1D::~CrankNicholson1D() { }
void CrankNicholson1D::preprt() {

}

void CrankNicholson1D::genframes(const vector<complex<double>>& v) {
    // first we should prepare Host and device vectors
    cuDoubleComplex* vecl = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex) * nx);
    cuDoubleComplex* vecr = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex) * nx);
    cuDoubleComplex* dev_vecl;
    cuDoubleComplex* dev_vecr;

    // Initialize the cuSparse context
    cusparseHandle_t h;
    cusparseStatus_t status= cusparseCreate(&h);
    // check for errors!
    if(status != CUSPARSE_STATUS_SUCCESS )
       DEBUG("Init failed!");

    //  after given an initial state, we initialize out internally used cuDoubleComplex vectors
    for(auto i = v.begin(); i != v.end(); ++i) {
        // the first vector is always on the right hand side first since
        // it's psi(t=0)
        auto j = distance( i, v.begin());
        vecr[j].x = i->real();
        vecr[j].y = i->imag();
    }

    // Initialize Cuda Memory/cusparse
    cudaMalloc((void**) dev_vecl, sizeof(cuDoubleComplex)*nx);
    cudaMalloc((void**) dev_vecr, sizeof(cuDoubleComplex)*nx);
    cudaMemcpy(dev_vecr,vecr,sizeof(cuDoubleComplex)*nx,cudaMemcpyHostToDevice);

    // Next on we initialize the 6 Arrays the the triangular Matrix
    // Only preperation of left hand side is necessary since the
    // right hand side can be calculated more efficiently
    cuDoubleComplex* opl1 = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex) * nx);
    cuDoubleComplex* opl2 = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex) * nx);
    cuDoubleComplex* opl3 = (cuDoubleComplex*) malloc(sizeof(cuDoubleComplex) * nx);
    cuDoubleComplex* opl1_dev;
    cuDoubleComplex* opl2_dev;
    cuDoubleComplex* opl3_dev;
    cudaMalloc((void**) &opl1_dev, sizeof(cuDoubleComplex) * nx);
    cudaMalloc((void**) &opl2_dev, sizeof(cuDoubleComplex) * nx);
    cudaMalloc((void**) &opl3_dev, sizeof(cuDoubleComplex) * nx);
    

}