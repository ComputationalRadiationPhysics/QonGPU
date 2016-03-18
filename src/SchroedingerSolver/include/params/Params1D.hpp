#ifndef PARAMS1D_H
#define PARAMS1D_H
using namespace std;

class Params1D: Params{

public:
	Params1D(complex<double> xma,complex<double> xmi,complex<double> tma,complex<double> tmi,int gridX,int gridT, int gridE): xmax(xma),xmin(xmi),tmax(tma),tmin(tmi),nx(gridX),nt(gridT),ne(gridE){};
private:
complex<double> xmax;
complex<double> xmin;
complex<double> tmax;
complex<double> tmin;
int nx,nt,ne;
};

#endif
