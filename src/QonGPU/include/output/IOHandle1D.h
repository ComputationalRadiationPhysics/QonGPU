//
// Created by zaph0d on 21/04/16.
//
#pragma once

#include "hdf5.h"
#include "hdf5_hl.h"
#include <thrust/host_vector.h>

class IOHandle1D {

public:

    IOHandle1D(Params1D* _p);
    ~IOHandle1D();
    cache_flush( const thrust::host_vector<cuDoubleComplex>& v);
private:

    char* filename = "framedata.h5";
    const size_t c_size = 100;
    thrust::host_vector<double> real_cache;
    thrust::host_vector<double> img_cache;
    // Require all the nececarry handles
    hid_t file;
    hid_t dataspace, dataset;
    hid_t filespace, memspace;
    hid_t prop;
    // Also size information will be
    // necessary
    hsize_t size;
    hsize_t offset;
    hsize_t dim1;


    void overwrite( const thrust::host_vector<cuDoubleComplex>& v);
    void savechunkt();
    void initfile();
};