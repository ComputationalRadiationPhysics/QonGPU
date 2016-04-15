#ifndef PARAMS1D_H
#define PARAMS1D_H
using namespace std;

class Params1D: Params{

public:
	Params1D(complex<double> xma,complex<double> xmi,complex<double> tma,complex<double> tmi,int gridX,int gridT, int gridE,double zin)
: xmax(xma),xmin(xmi),tmax(tma),tmin(tmi),nx(gridX),nt(gridT),ne(gridE),z(zin){};
    size_t getnx();
    size_t getne();
    size_t getnt();
    double getxmax();
    double getxmin();
    double gettmax();
    double gettmin();
    void setz(double zn);
    double getz();
private:
    complex<double> xmax;
    complex<double> xmin;
    complex<double> tmax;
    complex<double> tmin;
    size_t nx,nt,ne;
    complex<double> z;
};

double Params1D::getxmax(){
	return xmax.real();
}
double Params1D::getxmin(){return xmin.real();}
double Params1D::gettmax(){return tmax.real();}
double Params1D::gettmin(){return tmin.real();}
double Params1D::getz(){return z.real();}
size_t Params1D::getne(){return ne;}
size_t Params1D::getnx(){return nx;}
size_t Params1D::getnt(){return nt;}
void Params1D::setz(double zn){ z=zn;}
#endif
