//
// Created by zaph0d on 21/04/16.
//

#ifndef IOHANDLE1D_H
#define IOHANDLE1D_H

#include <vector>
#include "hdf5.h"
#include "hdf5_hl.h"
// We'll use  a struct of arrays to
// handle complex numbers in HDF5

#define FILENAME "simdata.h5"
#define DATASETNAME "Psi"
#define RANK 1
#define TCHUNK 1000
typedef struct comp{
    double *r;
    double *i;
    comp(hsize_t nx);
    comp(){};
    ~comp();
} comp;
comp::comp(hsize_t nx) {
    r = (double*) malloc(sizeof(double)*nx*TCHUNK);
    i = (double*) malloc(sizeof(double)*nx*TCHUNK);
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
    const size_t  nx,ne,nt;
    Params1D *param;
    void initfile();
    void closefile();
    void savechunk(const std::vector<std::complex<double>>& c);
    void copychunk(const std::vector<std::complex<double>>& c);

    // Some necessary HDF5 Variables
    hid_t file;
    hid_t dataset;
    hid_t memspace;
    hid_t filespace;
    hid_t prop;

    const hsize_t maxdims = H5F_UNLIMITED;
    hsize_t dset_dims = nx*ne;
    hsize_t chunk_dims = nx*TCHUNK;
};


#endif