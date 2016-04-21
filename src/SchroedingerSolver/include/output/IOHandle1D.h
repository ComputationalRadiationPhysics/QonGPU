//
// Created by zaph0d on 21/04/16.
//

#ifndef IOHANDLE1D_H
#define IOHANDLE1D_H

#define FILENAME "simdata.h5"
#define DATASETNAME "Psi"

#include "hdf5.h"
#include "hdf5_hl.h"
// We'll use  a struct of arrays to
// handle complex numbers in HDF5
typedef struct comp{
    double *r;
    double *i;
    comp(hsize_t nx);
    comp(){};
    ~comp();
} comp;
comp::comp(hsize_t nx) {
    r = (double*) malloc(sizeof(double)*nx);
    i = (double*) malloc(sizeof(double)*nx);
}
comp::~comp(){
    free(r);
    free(i);
}

class IOHandle1D {
public:
    IOHandle1D();
    IOHandle1D(Params1D *p);
    ~IOHandle1D();
private:
    comp chunk;
    hid_t file;
    const size_t  nx,ne,nt;
    Params1D *param;
    void initfile();
    void initdset();
    void closefile();
    void savechunk();
    
};


#endif