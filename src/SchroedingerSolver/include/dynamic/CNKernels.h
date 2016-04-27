//
// Created by zaph0d on 27/04/16.
//

#ifndef GIT_CNKERNELS_H
#define GIT_CNKERNELS_H

device_vector<cuDoubleComplex> operator+(device_vector<cuDoubleComplex> a, device_vector<cuDoubleComplex> b) {
    for(auto i = 0; i < a.size(); ++i) {
        a[i] += b[i];
    }
    return a;
}

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
    const double h = ( xmax - x0) / ( double) nx;
    const cuDoubleComplex h1 = make_cuDoubleComplex( -1.0/(h*h), 0);
    const cuDoubleComplex h2 = make_cuDoubleComplex( tau/2,0);
    double x = 0;

    while(ind < nx) {
        x = x0 + h*ind;
        rhs_rt( &in[ind], &in[ind-1], &in[ind+1], &out[ind], &s1, &s2, h1, h2, x);
        ind += oset;
    }
}

__global__ void create_const_diag(cuDoubleComplex* dl, cuDoubleComplex* du, double c, size_t nx) {
    // ensure that first element is not overwritten
    int tid = 1 +  threadIdx.x + blockDim.x * blockIdx.x;
    int oset = blockDim.x * gridDim.x;
    du[nx] = make_cuDoubleComplex(0,0);
    dl[nx] = make_cuDoubleComplex(0,0);
    cuDoubleComplex cc = make_cuDoubleComplex(0,c);
    while( tid < nx-1) {
        du[tid] = cc;
        dl[tid] = cc;
        tid += oset;
    }
}

__device__ __host__ inline void transform_diag( cuDoubleComplex *d, cuDoubleComplex* s1, cuDoubleComplex* s2,const double c, const double x, const cuDoubleComplex t1) {

    *s2 = make_cuDoubleComplex( -2 * c, 0) + pot(x);
    *s2 = *(s2) * t1;
    *s2 = make_cuDoubleComplex( -s2->y, s2->x);
    *d = *s1 + *s2;

}

__global__ void update_diagl( cuDoubleComplex* d, const double tau, const double h, const double c,const double xmin, const size_t nx, double t){

    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    int oset = blockDim.x * gridDim.x;
    auto t1 = make_cuDoubleComplex(tau/2 ,0);
    double x;
    cuDoubleComplex s1 = make_cuDoubleComplex(1.0 , 0);
    cuDoubleComplex s2;
    while( tid < nx) {
        x = xmin + h * (double) nx;
        transform_diag( &d[tid], &s1, &s2, c, x, t1);
        tid += oset;
    }

}

#endif
