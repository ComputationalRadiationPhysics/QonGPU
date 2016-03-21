#ifndef PARAMS1D_H
#define PARAMS1D_H
using namespace std;

class Params1D: Params{

public:
	Params1D(complex<double> xma,complex<double> xmi,complex<double> tma,complex<double> tmi,int gridX,int gridT, int gridE): xmax(xma),xmin(xmi),tmax(tma),tmin(tmi),nx(gridX),nt(gridT),ne(gridE){};
	int getnx();
	int getne();
	int getnt();
	double getxmax();
	double getxmin();
	double gettmax();
	double gettmin();
private:
complex<double> xmax;
complex<double> xmin;
complex<double> tmax;
complex<double> tmin;
int nx,nt,ne;

};

double Params1D::getxmax(){return xmax.real();}
double Params1D::getxmin(){return xmin.real();}
double Params1D::gettmax(){return tmax.real();}
double Params1D::gettmin(){return tmin.real();}
int Params1D::getne(){return ne;}
int Params1D::getnx(){return nx;}
int Params1D::getnt(){return nt;}
#endif
