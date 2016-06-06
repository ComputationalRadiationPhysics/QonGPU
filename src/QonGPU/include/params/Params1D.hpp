#pragma once

#include "Params.hpp"
using namespace std;

class Params1D: Params{

public:
	Params1D(double xma,
             double xmi,
             double tma,
             double tmi,
             int gridX,int gridT,
             int gridE,
             double zin,
             std::string name):
			 xmax(xma),
             xmin(xmi),
             tmax(tma),
             tmin(tmi),
             nx(gridX),
             nt(gridT),
             ne(gridE),
             z1(zin),
             filename(name){};

	size_t getnx();
    size_t getne();
    size_t getnt();
    double getxmax();
    double getxmin();
    double gettmax();
    double gettmin();
    void setz(double zn);
    double getz();
 	std::string getname() {return filename;};
	void seten(double _E) { E = _E;};
	double geten(){ return E;};
private:
	std::string filename;
    double xmax;
    double xmin;
    double tmax;
    double tmin;
    size_t nx,nt,ne;
    int x1,y1,z1;
    complex<double> z;
	double E  = 0;

};

double Params1D::getxmax(){
	return xmax;
}
double Params1D::getxmin(){ return xmin;}
double Params1D::gettmax(){ return tmax;}
double Params1D::gettmin(){ return tmin;}
double Params1D::getz(){ return z1;}
size_t Params1D::getne(){return ne;}
size_t Params1D::getnx(){return nx;}
size_t Params1D::getnt(){return nt;}
void Params1D::setz(double zn){ z=zn;}

