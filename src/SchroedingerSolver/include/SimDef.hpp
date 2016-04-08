#ifndef SRC_SCHROEDINGERSOLVER_INCLUDE_SIMDEF_HPP_
#define SRC_SCHROEDINGERSOLVER_INCLUDE_SIMDEF_HPP_

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
  explicit SimDef(Params1D *p):da(p), s(p, psi0) {
    DEBUG("CALL CONSTRUCTOR")
        }
  explicit SimDef(Params2D *p):da(p) {
  }
  explicit SimDef(Params3D *p):da(p) {
  }
  void staticsolve() {
    s.solve();
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
  std::complex<double>* psi0;
  std::complex<double>* psi;
  double* npsi;
};


#endif
