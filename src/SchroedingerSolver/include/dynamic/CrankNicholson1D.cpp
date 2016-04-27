//
// Created by max on 22/04/16.
//

#include "CrankNicholson1D.hpp"

__device__ __host__ cuDoubleComplex inline pot(double x){
    return make_cuDoubleComplex(1/sqrt(x*x+1),0);
}

__device__ __host__ inline void rhs_rt(
                                         cuDoubleComplex* in1,
                                         cuDoubleComplex* in2,
                                         cuDoubleComplex* in3,
                                         cuDoubleComplex* out,
                                         cuDoubleComplex* s1,
                                         cuDoubleComplex* s2,
                                         const cuDoubleComplex& h1,
                                         const cuDoubleComplex& h2,
                                         double x
                                         ) {


    *s1 = h1 * ( *in2 + *(in3)  - make_cuDoubleComplex( 2.0, 0) * *in1);
    *s2 = pot(x) * *in1;
    *s1 = *s1 + *s2;
    *s1 = make_cuDoubleComplex( -s1->y, s2->x);
    *out = *s1 + *in1;
}

__global__ void transform_rhs(cuDoubleComplex* in, cuDoubleComplex* out,size_t nx,double xmax,double x0,double tau) {

    int ind = threadIdx.x + blockDim.x * blockIdx.x;
    int oset = blockDim.x * gridDim.x;
    cuDoubleComplex s1,s2;
    const cuDoubleComplex h1 = make_cuDoubleComplex( -1.0, 0);
    const cuDoubleComplex h2 = make_cuDoubleComplex( tau/2,0);
    const double h = ( xmax - x0) / ( double) nx;
    double x = 0;
    while(ind < nx) {
        // Those lines need to be tested first!
        x = x0 + h*ind;
        rhs_rt( &in[ind], &in[ind-1], &in[ind+1], &out[ind], &s1, &s2, h1, h2, x);
        ind += oset;
        /* x1 = rhs[ind];
        x2 = rhs[ind-1];
        x3 = rhs[ind+1];
        s1 = h1*(x2 + x3 - make_cuDoubleComplex(2,0) * x1);
        s2 = pot(x)*x1;
        s1 = s1 + s2;
        s1 = make_cuDoubleComplex(-s1.y,s2.x);
        rhs[ind] = rhs[ind] + s1;
        ind += oset; */
    }
}

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
                                                                               xmin( _p->getxmin())
{
    thrust::copy( chunk_h.begin(), chunk_h.end(), chunkr_d.begin());
}
CrankNicholson1D::~CrankNicholson1D() { }

void CrankNicholson1D::rhs_rt( double x, double tau) {

    transform_rhs(raw_pointer_cast(chunkl_d.data()), raw_pointer_cast(chunkr_d.data()), nx, param->getxmax(), param->getxmin(),tau);
}

void CrankNicholson1D::lhs_rt(double x, double t,
                              cuDoubleComplex* d,
                              cuDoubleComplex* ud,
                              cuDoubleComplex* du ) {


}