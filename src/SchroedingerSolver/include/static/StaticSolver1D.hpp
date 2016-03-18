#ifndef STATICSOLVER1D_H
#define STATICSOLVER1D_H

class StaticSolver1D: StaticSolver {
public:
	StaticSolver1D(Params1D& p,std::complex<double>* ps):psi0(ps){};
	StaticSolver1D(){};
    void solve();
protected:
	std::complex<double>* psi0;
    Params1D* p;
};


#endif 
