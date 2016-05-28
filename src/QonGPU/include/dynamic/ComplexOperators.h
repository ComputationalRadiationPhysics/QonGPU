/*
 * This class has to be seen as an adapter
 * to use cuDoublecComplex with Operators in a
 * more natural way!
 *
 */
#ifndef COMPLEXOPERATORS_H
#define COMPLEXOPERATORS_H

#include "cuComplex.h"

struct add {
   __host__ __device__ __inline__ cuDoubleComplex operator()(cuDoubleComplex s1, cuDoubleComplex s2) {
        return cuCadd(s1, s2);
    }
};

struct  sub {
    __host__ __device__ __inline__  cuDoubleComplex operator()(cuDoubleComplex m1, cuDoubleComplex m2) {
        return cuCsub(m1,m2);
    }
};

struct mult {
    __host__ __device__ __inline__ cuDoubleComplex operator()(cuDoubleComplex p1, cuDoubleComplex p2) {
        return cuCmul(p1,p2);
    }
};

struct div {
    __host__ __device__ __inline__ cuDoubleComplex operator()(cuDoubleComplex div1, cuDoubleComplex div2) {
        return cuCdiv(div1, div2);
    }
};

struct conj {
    __host__ __device__ __inline__ cuDoubleComplex operator()(cuDoubleComplex val) {
        return make_cuDoubleComplex(val.x, -val.y);
    }
};

struct  abs {
    __host__ __device__ __inline__ double operator()(cuDoubleComplex val) {
        return cuCabs(val);
    }
};


__host__ __device__ __inline__ cuDoubleComplex operator+(cuDoubleComplex a,const cuDoubleComplex& b) {
    return cuCadd(a,b);
}


__host__ __device__ __inline__ cuDoubleComplex operator-(cuDoubleComplex a, const cuDoubleComplex& b) {
    return cuCsub(a,b);
}


__host__ __device__ __inline__ cuDoubleComplex operator*(cuDoubleComplex a, cuDoubleComplex b) {
    return cuCmul(a,b);
}


__host__ __device__ __inline__ cuDoubleComplex operator/(cuDoubleComplex a, cuDoubleComplex b) {
    return cuCdiv(a,b);
}


__host__ __device__ __inline__ cuDoubleComplex operator+=(cuDoubleComplex a, cuDoubleComplex b){
    return a+b;
}
__host__ __device__ __inline__ cuDoubleComplex operator*(const double a, cuDoubleComplex b) {

    return make_cuDoubleComplex(b.x*a, a*b.y);
}


std::ostream& operator<<( std::ostream &out, const cuDoubleComplex z) {

    out<<" "<< z.x << " " << z.y;

    return out;
}
#endif