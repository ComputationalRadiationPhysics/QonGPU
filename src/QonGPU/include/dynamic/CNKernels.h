//
// Created by zaph0d on 27/04/16.
//
#pragma once

#ifndef CUDART_PI_F
#define CUDART_PI_F 3.141592653589793
#endif

__device__ __host__ cuDoubleComplex pot(double x, double t) {
	
	
	double a = 1.36917961592e-07;
	double b = 0.0069314718056;
	double t0 = 2500.0;
	double w = 0.08607963870836033;
	double k = w/137;
	double I = 0.3;
	// Only have time-dependence if t>0
	
	
	double g1 = exp(-a*(t-t0)*(t-t0));
	double g2 = exp(-b*x*x);
	double f = pow(sin(w*t - k*x),1);
	double res = -1/sqrt(x*x+1) + f*I*g1*g2;
	
	return make_cuDoubleComplex(res, 0);
    
}


__device__ __host__ inline void mult_rhs( cuDoubleComplex& in1,
                                          cuDoubleComplex& in2,
                                          cuDoubleComplex& in3,
                                          cuDoubleComplex& out,
                                          const cuDoubleComplex& h1,
                                          const cuDoubleComplex& h2,
                                          double x,
                                          double t) {


    //s1 = (h1 * h2) * ( (in2) + (in3)  + make_cuDoubleComplex(-2 * in1.x, 2 * in1.y));

    cuDoubleComplex s1 = make_cuDoubleComplex( 0, 0);
    cuDoubleComplex s2 = make_cuDoubleComplex( 0, 0);

    s1 = cuCadd(in2 , in3);
    s1 = cuCadd(s1, make_cuDoubleComplex(-2 * in1.x, -2 * in1.y));
    s1 = cuCmul(h1, s1);
    s1 = cuCmul(h2, s1);

    s2 = cuCmul(h2, pot(x, t));
    s2 = cuCmul(s2, in1);

    s1 = cuCadd(s1, s2);
    s1 = make_cuDoubleComplex( -s1.y, s1.x);
    out =  cuCadd(in1, s1) ;

}


__global__ void transform_rhs(cuDoubleComplex* in, // note that in is just an temporary array to
                              cuDoubleComplex* out,// ensure threat safety!
                              size_t nx,
                              double xmax,
                              double xmin,
                              double tau,
                              double t) {

    // This function multiplies the psi(t)
    // with the Crank Nicholson time
    // Time evolution operator
    int ind = threadIdx.x + blockDim.x * blockIdx.x;
    int oset = blockDim.x * gridDim.x;
    const double h = ( xmax - xmin) / ( double) nx;
    const cuDoubleComplex h1 = make_cuDoubleComplex( -1.0 / ( 2.0 * h * h), 0);
    const cuDoubleComplex h2 = make_cuDoubleComplex( -tau / 2.0 ,0);
    double x = xmin;

    if(ind == 0 ) {

        cuDoubleComplex bound = make_cuDoubleComplex(0.0, 0.0);
        mult_rhs(in[ind], bound, in[ind + 1], out[ind], h1, h2, x, t);
        ind += oset;

    }

        while (ind < nx) {
            x = xmin + h * (double) ind;
            if(ind == nx-1) {

                cuDoubleComplex s1 = make_cuDoubleComplex( 0, 0);
                cuDoubleComplex s2 = make_cuDoubleComplex( 0, 0);

                cuDoubleComplex  bound = make_cuDoubleComplex(0,0);
                s1 = in[nx-2];
                s1 = cuCadd(s1, make_cuDoubleComplex(-2* in[ind].x, -2 * in[ind].y));
                s1 = cuCmul(h1,s1);
                s1= cuCmul(h2, s1);

                s2 = cuCmul(h2, pot(x, t));
                s2 = cuCmul(s2, in[ind]);
                s1 = cuCadd(s1,s2);
                s1 = make_cuDoubleComplex( -s1.y, s1.x);
                out[ind] = cuCadd(in[ind], s1);
                ind += oset;

            }
            else {
                mult_rhs(in[ind], in[ind - 1], in[ind + 1], out[ind],
                          h1, h2, x, t);

            }
            //printf("out = %lf + i %lf \n", out[ind].x, out[ind].y);
            ind += oset;

        }


}


__global__ void create_const_diag(	cuDoubleComplex* dl,
									cuDoubleComplex* du,
									double c,
									size_t nx ) {


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
                                                const cuDoubleComplex& t1,
                                                double t) {

    cuDoubleComplex temp1 = make_cuDoubleComplex( c, 0);
    s2 = cuCadd(temp1, pot(x, t));
    s2 = cuCmul(s2, t1);
    s2 = make_cuDoubleComplex( - s2.y , s2.x);
    d = cuCadd(s2, s1);

}


__global__ void update_diagl( cuDoubleComplex* d,
                              const double tau,
                              const double h,
                              const double xmin,
                              const size_t nx,
                              double t) {

    int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid == 1000)
			printf("Diagonal update called! \n");

    int oset = blockDim.x * gridDim.x;
    cuDoubleComplex t1 = make_cuDoubleComplex( tau / 2.0 ,0);

    // 2  and  - is left out since, you can make the
    // expression easier by that!
    double c =  1.0 / (h * h);
    double x = xmin;
    cuDoubleComplex s1 = make_cuDoubleComplex( 1.0 , 0);
    cuDoubleComplex s2 = make_cuDoubleComplex( 0, 0);

    while( tid < nx) {

        x = xmin + h * (double) tid;
        transform_diag( d[tid], s1, s2, c, x, t1, t);
        tid += oset;

    }

}
