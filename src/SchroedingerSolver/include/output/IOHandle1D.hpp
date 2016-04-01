#ifndef IOHANDLE1D_H_
#define IOHANDLE1D_H_

#define FILENAME "elevels.h5"

#include "IOHandle.hpp"

#include <vector>

#include <string>
using namespace std;
class IOHandle1D:IOHandle {
public:
    IOHandle1D(){};
    IOHandle1D(Params1D*);
    ~IOHandle1D();
    void savelevel(vector<double> psi,vector<double> e,int z);
private:
    Params1D* p;
    hid_t file;
};

IOHandle1D::IOHandle1D(Params1D* pa):p(pa) {
    file=H5Fcreate(FILENAME,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    double* tempdata = (double*) malloc(sizeof(double)*7);
    hsize_t dim=7;
    tempdata[0]=p->getxmax();
    tempdata[1]=p->getxmin();
    tempdata[2]=p->gettmax();
    tempdata[3]=p->gettmin();
    tempdata[4]=(double) p->getnx();
    tempdata[5]=(double) p->getnt();
    tempdata[6]=(double) p->getne();
    H5LTmake_dataset(file,"/params1d",1,&dim, H5T_NATIVE_DOUBLE,tempdata);
    free(tempdata);
}
IOHandle1D::~IOHandle1D() {
    H5Fclose(file);
}
void IOHandle1D::savelevel(vector<double> psi,vector<double> e,int z) {
    
}
#endif
