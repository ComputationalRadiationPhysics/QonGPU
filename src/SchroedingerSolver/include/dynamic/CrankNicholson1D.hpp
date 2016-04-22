//
// Created by max on 20/04/16.
//


#ifndef PROJECT_CRANKNICHOLSON1D_HPP
#define PROJECT_CRANKNICHOLSON1D_HPP

#include <vector>
#include <complex>
#include <cuda_runtime.h>
#include <cusparse_v2.h>

#include "TimeOperator1D.hpp"

#define TCHUNK 100







using namespace std;
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
    CrankNicholson1D(Params1D& p);
    ~CrankNicholson1D();
    void genframes(const vector<complex<double>>& v);
private:
    const size_t nx,nt;
    Params1D* params;
    // Cuda uses PODs so we'll internally stick with pods!
    vector<complex<double>> frames;
    void saveChunk();
    // Defines a preperation routine, memory routine, solve rountine
    void preprt(cuDoubleComplex* vecr_dev, cuDoubleComplex* opl1_dev, cuDoubleComplex* opl2_dev, cuDoubleComplex* opl3_dev);
    void memrt();
    void solvert();
};
#endif
