//
// Created by zaph0d on 27/04/16.
//
#pragma once

__device__ __host__ cuDoubleComplex pot(double x) {

    //return make_cuDoubleComplex( - 1.0 /sqrt(x * x + 1.0), 0);
    return make_cuDoubleComplex(x*x/2.0, 0);
}


__device__ __host__ inline void mult_rhs( cuDoubleComplex& in1,
                                          cuDoubleComplex& in2,
                                          cuDoubleComplex& in3,
                                          cuDoubleComplex& out,
                                          const cuDoubleComplex& h1,
                                          const cuDoubleComplex& h2,
                                          double x) {


    //s1 = (h1 * h2) * ( (in2) + (in3)  - make_cuDoubleComplex(2 * in1.x, 2 * in1.y));

    cuDoubleComplex s1 = make_cuDoubleComplex( 0, 0);
    cuDoubleComplex s2 = make_cuDoubleComplex( 0, 0);

    s1 = cuCadd(in2 , in3);
    s1 = cuCadd(s1, make_cuDoubleComplex(-2 * in1.x, -2 * in1.y));
    s1 = cuCmul(h1, s1);
    s1 = cuCmul(h2, s1);

    s2 = cuCmul(h2, pot(x));
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
                              double tau) {

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
        mult_rhs(in[ind], bound, in[ind + 1], out[ind], h1, h2, x);
        ind += oset;

    }

    if(ind == 1) {
        x = xmin + h;
        //mult_rhs(in[1], in[0], in[2], out[1], h1, h2, x);
        cuDoubleComplex s1 = make_cuDoubleComplex( 0, 0);
        cuDoubleComplex s2 = make_cuDoubleComplex( 0, 0);

        s1 = cuCadd(in[0] , in[2]);
        //printf("s1 = %f + i %lf \n", s1.x, s1.y);
        s1 = cuCadd(s1, make_cuDoubleComplex(-2 * in[1].x, -2 * in[1].y));
        //printf("s1 = %f + i %lf \n", s1.x, s1.y);
        s1 = cuCmul(h1, s1);
        //printf("s1 = %f + i %lf \n", s1.x, s1.y);
        s1 = cuCmul(h2, s1);
        //printf("s1 = %f + i %lf \n", s1.x, s1.y);


        s2 = cuCmul(h2, pot(x));
        //printf("s2 = %f + i %lf \n", s2.x, s2.y);
        s2 = cuCmul(s2, in[1]);
        //printf("s2 = %f + i %lf \n", s2.x, s2.y);

        s1 = cuCadd(s1, s2);
        //printf("s1 = %f + i %lf \n", s1.x, s1.y);

        s1 = make_cuDoubleComplex( -s1.y, s1.x);
        //printf("s1 = %f + i %lf \n", s1.x, s1.y);
        out[1] =  cuCadd(in[1], s1);
        //printf("out = %f + i %lf \n", out[1].x, out[1].y);

        ind +=oset;
        //printf("h1 = %f + i %f \n", h1.x, h1.y);
        //printf("h2 = %f + i %f \n", h2.x, h2.y);
        //printf("The result: out = %f + i %f \n", out[1].x, out[1].y);
        //printf("in[2] = %f + i %f \n", in[2].x, in[2].y);
        //printf("in[1] = %f + i %f \n", in[1].x, in[1].y);
        //printf("in[0] = %f + i %f \n", in[0].x, in[0].y);
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

                s2 = cuCmul(h2, pot(x));
                s2 = cuCmul(s2, in[ind]);
                s1 = cuCadd(s1,s2);
                s1 = make_cuDoubleComplex( -s1.y, s1.x);
                out[ind] = cuCadd(in[ind], s1);
                ind += oset;

            }
            else {
                mult_rhs(in[ind], in[ind - 1], in[ind + 1], out[ind],
                          h1, h2, x);

            }
            //printf("out = %lf + i %lf \n", out[ind].x, out[ind].y);
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

    cuDoubleComplex temp1 = make_cuDoubleComplex( c, 0);
    s2 = cuCadd(temp1, pot(x));
    s2 = cuCmul(s2, t1);
    s2 = make_cuDoubleComplex( - s2.y , s2.x);
    d = cuCadd(s2, s1);
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
    cuDoubleComplex t1 = make_cuDoubleComplex( tau / 2.0 ,0);

    // 2  and  - is left out since, you can make the
    // expression easier by that!
    double c =  1.0 / (h * h);
    double x = xmin;
    cuDoubleComplex s1 = make_cuDoubleComplex( 1.0 , 0);
    cuDoubleComplex s2 = make_cuDoubleComplex( 0, 0);

    while( tid < nx) {

        x = xmin + h * (double) tid;
        transform_diag( d[tid], s1, s2, c, x, t1);
        tid += oset;

    }

}
