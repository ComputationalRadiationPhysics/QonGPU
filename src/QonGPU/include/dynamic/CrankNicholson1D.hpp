#pragma once

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <mpi.h>


#include "splash/splash.h"
#include "../output/IOHandle1D.cpp"
#include "TimeOperator1D.hpp"
#include "ComplexOperators.h"
#include "spike_host.cu"

using namespace thrust;
class CrankNicholson1D : TimeOperator1D {
    /** This class implements the Crank Nicholson
     *  algorithm for the 1Dim Schroedinger Equation
     *  in Cuda C. The idea is that each Timestep
     *  calculation is done by using the cusparse
     *  Tridiagonal solver Algorithm(which is hopefully
     *  fast/accurate enough for our application)
     */
public:
    CrankNicholson1D(Params1D *_p);
    ~CrankNicholson1D();
    void time_solve();
    size_t getnx(){ return nx;};
    size_t getnt(){ return nt;}
    double gettmax(){ return tmax;};
    double gettmin(){ return tmin;};
    double getxmax(){ return xmax;};
    double getxmin(){ return xmin;};
    void cusparse_init();
    void cusparse_destr();
    // get a function to copy the initial state!
    void setstate(const thrust::host_vector<cuDoubleComplex>& v);

private:
    Params1D* param;
    IOHandle1D io;
    const size_t nx,nt;
    const double E;
    const double tau;
    const size_t csize = 100;
    double tmax, xmax;
    double tmin, xmin;
    host_vector<cuDoubleComplex> inital;
    host_vector<cuDoubleComplex> chunk_h;
    // lefthand side chunk on Device
    device_vector<cuDoubleComplex>  chunkl_d;
    // righthand side chunk on Device
    device_vector<cuDoubleComplex> chunkr_d;
    // Diagonals of the triangular matrix
    device_vector<cuDoubleComplex> d, dl, du;
    splash::SerialDataCollector HDFile;



    std::string filename;

    // Define necessary members
    cusparseStatus_t status, status2;
    cusparseHandle_t handle = 0;
    // Define cusparse member functions
    // @TODO export this stuff into an interface

    void cusparse_sv();

    // Define necessary member functions
    void rhs_rt( const double c);
    void lhs_rt( double x, double t,
                 cuDoubleComplex* d,
                 cuDoubleComplex* du,
                 cuDoubleComplex* dl);
    void prp_rt();
    void printinitial();
    void initfile(splash::DataCollector::FileCreationAttr& fa);
    void savechunk(int step);
    void closefile();
    void save_vector(int step, const thrust::device_vector<cuDoubleComplex>& v);
    void save_vectorh(int step, const thrust::host_vector<cuDoubleComplex>& v);
    void save_diag(int step, const thrust::device_vector<cuDoubleComplex>& v);


};



