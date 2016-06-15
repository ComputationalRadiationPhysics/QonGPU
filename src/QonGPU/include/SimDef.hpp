#pragma  once
#include <cuda_runtime.h>
#include <cublas_v2.h>

#include "cuComplex.h"


#include "AllHeader.hpp"


template <class StatSolver, class TimeOp, class Pot,
          class Para, int dim>
class SimDef{
  
public:

    SimDef(Params1D *p):da(p), s(p),t(p), psi0(p->getnx()){
       
    }

    SimDef(Params2D *p):da(p) {}
    SimDef(Params3D *p):da(p) {}

    void staticsolve() {
        s.solve();
        s.copystate(0, psi0);
        t.setstate(psi0);
       
    }

    void timerev() {

        std::chrono::high_resolution_clock::time_point t1 =
                std::chrono::high_resolution_clock::now();
        DEBUG2("Time Solve");
        t.time_solve();
        std::chrono::high_resolution_clock::time_point t2 =
                std::chrono::high_resolution_clock::now();
        auto duration =
                std::chrono::duration_cast<std::chrono::minutes>( t2 - t1 ).count();

        std::cout<<"It took "<< duration <<" minutes for the whole simulation!"<<std::endl;

    }


private:

    StatSolver s;
    TimeOp t;
    Pot p;
    Para *da;
    thrust::host_vector<cuDoubleComplex> psi0;
};
