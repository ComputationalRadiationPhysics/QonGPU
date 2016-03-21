#ifndef STATICSOLVER1D_H
#define STATICSOLVER1D_H

class StaticSolver1D: StaticSolver {
public:
	StaticSolver1D(Params1D* pa,std::complex<double>* ps):p(pa),psi0(ps){};
	StaticSolver1D(){};
    void staticsolve();
protected:
	std::complex<double>* psi0;
    Params1D* p;
};


#endif 
