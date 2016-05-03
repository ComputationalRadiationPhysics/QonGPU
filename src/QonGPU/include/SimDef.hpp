#pragma  once
#include <cuda_runtime.h>

#include <cublas_v2.h>

#include "cuComplex.h"

#include "boost/static_assert.hpp"

#include "boost/type_traits/is_base_of.hpp"

#include "thrust/complex.h"

#include "AllHeader.hpp"



template <class StatSolver, class TimeOp, class IO, class Pot,
          class Para, int dim>
class SimDef{
  
 public:
  explicit SimDef(Params1D *p):da(p), s(p) {}
  explicit SimDef(Params2D *p):da(p) {}
  explicit SimDef(Params3D *p):da(p) {}
  void staticsolve() {
    s.solve();
    //s.copystate(0, &psi0);
    //t.setstatic(&psi0);
    }
  void timerev() {
  }
  void printres() {
  }

	
 private:
  StatSolver s;
  TimeOp t;
  IO io;
  Pot p;
  Para *da;
  thrust::host_vector<cuDoubleComplex> psi0;
};
