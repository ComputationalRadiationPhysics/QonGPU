//
// Created by zaph0d on 27/04/16.
//
#pragma once

device_vector<cuDoubleComplex> operator+(device_vector<cuDoubleComplex> a,
                                         device_vector<cuDoubleComplex> b) {

    for(int i = 0; i < a.size(); ++i) {

        a[i] += b[i];

    }

    return a;
}


__device__ __host__ cuDoubleComplex inline pot(double x){

    return make_cuDoubleComplex( - 2.0 /sqrt(x * x + 1.0), 0);

}


__device__ __host__ inline void mult_rhs( cuDoubleComplex* in1,
                                          cuDoubleComplex* in2,
                                          cuDoubleComplex* in3,
                                          cuDoubleComplex* out,
                                          cuDoubleComplex* s1,
                                          cuDoubleComplex* s2,
                                          const cuDoubleComplex& h1,
                                          const cuDoubleComplex& h2,
                                          double x) {

    *s1 = h1 * h2 * ( (*in2) + (*in3)  - make_cuDoubleComplex( 2.0, 0) * (*in1));
    *s2 = h2 * pot(x) * (*in1);
    *s1 = (*s1) + (*s2);
    *s1 = make_cuDoubleComplex( -s1->y, s1->x);
    *out =  (*in1) + (*s1) ;

}


__global__ void transform_rhs(cuDoubleComplex* in, // note that in is just an temporary array to
                              cuDoubleComplex* out,// ensure threat safety!
                              size_t nx,
                              double xmax,
                              double xmin,
                              double tau, double c) {

    // This function multiplies the psi(t)
    // with the Crank Nicholson time
    // Time evolution operator
    int ind = threadIdx.x + blockDim.x * blockIdx.x;
    int oset = blockDim.x * gridDim.x;
    cuDoubleComplex s1 = make_cuDoubleComplex( 0, 0);
    cuDoubleComplex s2 = make_cuDoubleComplex( 0, 0);
    const double h = ( xmax - xmin) / ( double) nx;
    const cuDoubleComplex h1 = make_cuDoubleComplex( -1.0 / ( 2.0 * h * h), 0);
    const cuDoubleComplex h2 = make_cuDoubleComplex( -tau / 2.0 ,0);
    double x = xmin;

    while(ind < nx) {

        x +=  h * (double) ind;
        mult_rhs( &in[ind], &in[ind-1], &in[ind+1], &out[ind],
                  &s1, &s2, h1, h2, x);
        ind += oset;
    }

}


__global__ void create_const_diag(cuDoubleComplex* dl,
                                  cuDoubleComplex* du,
                                  double c,
                                  size_t nx) {


    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    int oset = blockDim.x * gridDim.x;
    cuDoubleComplex cc = make_cuDoubleComplex(0,c);

    while( tid < nx) {

        du[tid] = cc;
        dl[tid] = cc;
        tid += oset;

    }

    du[nx-1] = make_cuDoubleComplex(0,0);
    dl[0] = make_cuDoubleComplex(0,0);

}


__device__ __host__ inline void transform_diag( cuDoubleComplex& d,
                                                cuDoubleComplex& s1,
                                                cuDoubleComplex& s2,
                                                const double c,
                                                const double x,
                                                const cuDoubleComplex& t1) {

    s2 = make_cuDoubleComplex( c , 0) + pot(x);
    s2 = s2 * t1;
    s2 = make_cuDoubleComplex( s2.y * (-1.0), s2.x);
    d = s1 + s2;
    //printf("pot(x) = %lf \n", pot(x).x);
    //printf("D = %lf \n", d.y);

}


__global__ void update_diagl( cuDoubleComplex* d,
                              const double tau,
                              const double h,
                              const double xmin,
                              const size_t nx) {

    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    int oset = blockDim.x * gridDim.x;
    auto t1 = make_cuDoubleComplex( tau / 2.0 ,0);

    // 2  and  - is left out since, you can make the
    // expression easier by that!
    double c =  1.0 / ( h * h);
    double x = xmin;
    cuDoubleComplex s1 = make_cuDoubleComplex( 1.0 , 0);
    cuDoubleComplex s2 = make_cuDoubleComplex( 0, 0);

    while( tid < nx) {

        x += h * (double) tid;
        transform_diag( d[tid], s1, s2, c, x, t1);
        tid += oset;

    }

}
