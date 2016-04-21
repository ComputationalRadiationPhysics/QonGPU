//
// Created by max on 20/04/16.
//


#ifndef PROJECT_CRANKNICHOLSON1D_HPP
#define PROJECT_CRANKNICHOLSON1D_HPP



#include <vector>
#include <complex>

#include "TimeOperator1D.hpp"

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
private:
    const size_t nx,nt;
    Params1D* params;
    vector<complex<double>> psi0;

};
#endif
