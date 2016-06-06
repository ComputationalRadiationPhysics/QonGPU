#pragma once

#ifndef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#endif
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#ifndef __GCC__
#pragma GCC diagnostic pop
#endif

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <chrono>
#include <cuda_runtime.h>
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif
#include <cmath>

#include "ComplexOperators.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "cusparse.h"
#include "../output/IOHandle1D.cpp"
#include "TimeOperator1D.hpp"
#include "spike_host.cu"
#include "ThomasSerial.h"
#include "CNKernels.h"
#include "MatrixGeneration.h"
#include "TridiagMult.h"



using namespace thrust;

class CrankNicolson1D : TimeOperator1D {

    /** This class implements the Crank Nicholson
     *  algorithm for the 1Dim Schroedinger Equation
     *  in Cuda C. The idea is that each Timestep
     *  calculation is done by using the cusparse/SpikeCR
     *  Tridiagonal solver Algorithm(which is hopefully
     *  fast/accurate enough for our application)
     */

public:

    // Define constructor und destructor
    CrankNicolson1D(Params1D *_p);
    ~CrankNicolson1D();

    // Define Solver Function
    void time_solve();

    // get a function to copy the initial state!
    void setstate(const thrust::host_vector<cuDoubleComplex>& v);

private:

    Params1D* param;
    const size_t nx,nt;
    const double E;
    const double tau;
    const size_t csize = 100;
    double tmax, xmax;
    double tmin, xmin;


    // lefthand side chunk on Device
    device_vector<cuDoubleComplex>  chunkl_d;

    // righthand side chunk on Device
    device_vector<cuDoubleComplex> chunkr_d;

    // Diagonals of the triangular matrix
    device_vector<cuDoubleComplex> d, dl, du;


    // Filename from Command line input
    std::string filename;

    // Define necessary member functions
    void rhs_rt();
    void printinitial();
    void write_p(hid_t *f);

    // Define saving methods
    void savechunk(int step);
    void save_vector(int step, const thrust::device_vector<cuDoubleComplex>& v);
    void save_vectorh(int step, const thrust::host_vector<cuDoubleComplex>& v);
    void save_diag(int step, const thrust::device_vector<cuDoubleComplex>& v);


};



