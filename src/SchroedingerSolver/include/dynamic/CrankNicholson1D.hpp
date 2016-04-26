//
// Created by max on 20/04/16.
//


#ifndef PROJECT_CRANKNICHOLSON1D_HPP
#define PROJECT_CRANKNICHOLSON1D_HPP


#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cuda_runtime.h>
#include <cusparse_v2.h>

#include "TimeOperator1D.hpp"
#include "ComplexOperators.h"







using namespace thrust;
class CrankNicholson1D : TimeOperator1D {
    /** This class implements the Crank Nicholson
     *  algorithm for the 1Dim Schroedinger Equation
     *  in Cuda C. The idea is that each Timestep
     *  calculation is done by using the cusparest
     *  Triangular solver Algorithm(which is hopefully
     *  fast enough for our applications.
     */
public:
    CrankNicholson1D();
    CrankNicholson1D(Params1D *_p, vector<cuDoubleComplex> _v);
    ~CrankNicholson1D();
    void time_ev();
private:
    Params1D* param;
    const size_t nx,nt;
    const double E;
    const size_t csize = 100;
    host_vector<cuDoubleComplex> chunk_h;
    device_vector<cuDoubleComplex>  chunk_d;
    void rhs_rt();
    void lhs_rt();
    void prp_rt();
    void save_chunch();
};


device_vector<cuDoubleComplex> operator+(device_vector<cuDoubleComplex> a, device_vector<cuDoubleComplex> b) {
    for(auto i = 0; i < a.size(); ++i) {
        a[i] += b[i];
    }
    return a;
}

#endif
