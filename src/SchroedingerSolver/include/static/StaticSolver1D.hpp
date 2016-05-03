#ifndef STATICSOLVER1D_H
#define STATICSOLVER1D_H

class StaticSolver1D: StaticSolver {
public:
	StaticSolver1D(Params1D* pa):p(pa){};
	StaticSolver1D(){};
    void staticsolve();
protected:
    Params1D* p;
};


#endif 
