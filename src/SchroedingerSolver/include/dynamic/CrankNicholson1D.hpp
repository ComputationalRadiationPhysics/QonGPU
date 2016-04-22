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

#define CLEANUP(s)                                   \
do {                                                 \
    printf ("%s\n", s);                              \
    if (yHostPtr)           free(yHostPtr);          \
    if (zHostPtr)           free(zHostPtr);          \
    if (xIndHostPtr)        free(xIndHostPtr);       \
    if (xValHostPtr)        free(xValHostPtr);       \
    if (cooRowIndexHostPtr) free(cooRowIndexHostPtr);\
    if (cooColIndexHostPtr) free(cooColIndexHostPtr);\
    if (cooValHostPtr)      free(cooValHostPtr);     \
    if (y)                  cudaFree(y);             \
    if (z)                  cudaFree(z);             \
    if (xInd)               cudaFree(xInd);          \
    if (xVal)               cudaFree(xVal);          \
    if (csrRowPtr)          cudaFree(csrRowPtr);     \
    if (cooRowIndex)        cudaFree(cooRowIndex);   \
    if (cooColIndex)        cudaFree(cooColIndex);   \
    if (cooVal)             cudaFree(cooVal);        \
    if (descr)              cusparseDestroyMatDescr(descr);\
    if (handle)             cusparseDestroy(handle); \
    cudaDeviceReset();          \
    fflush (stdout);                                 \
} while (0)







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
    void preprt();
    void memrt();
    void solvert();
};
#endif
