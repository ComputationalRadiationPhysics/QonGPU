#ifndef SRC_SCHROEDINGERSOLVER_INCLUDE_STATIC_STATICSOLVER_H_
#define SRC_SCHROEDINGERSOLVER_INCLUDE_STATIC_STATICSOLVER_H_

class StaticSolver{
 public:
  StaticSolver(Params* p, double*  psi0) {}
  StaticSolver() {}
  void solve();
};

#endif




