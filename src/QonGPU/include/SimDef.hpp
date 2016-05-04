#pragma  once
#include <cuda_runtime.h>

#include <cublas_v2.h>

#include "cuComplex.h"

#include "thrust/complex.h"

#include "AllHeader.hpp"


template <class StatSolver, class TimeOp, class Pot,
          class Para, int dim>
class SimDef{
  
public:

    explicit SimDef(Params1D *p):da(p), s(p),t(p){

    }
    explicit SimDef(Params2D *p):da(p) {}
    explicit SimDef(Params3D *p):da(p) {}
    void staticsolve() {
        s.solve();
        s.copystate( 0,&psi0);
        t.setstate(psi0);
    }
    void timerev() {
        t.time_solve();
    }
    void printres() {
    }
private:

    StatSolver s;
    TimeOp t;
    Pot p;
    Para *da;
    thrust::host_vector<cuDoubleComplex> psi0;
};
