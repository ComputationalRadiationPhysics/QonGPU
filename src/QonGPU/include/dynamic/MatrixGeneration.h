//
// Created by zaph0d on 28/05/16.
//
#pragma once

#include "MemFunctors.h"
/*
void create_diags(thrust::device_vector<cuDoubleComplex>& dl,
                  thrust::device_vector<cuDoubleComplex>& du,
                  cuDoubleComplex c) {

    thrust::host_vector<cuDoubleComplex> vec( dl.size());
    fetch_vec f();
    push_vec p();


    //f(dl, vec);
    thrust::copy(dl.begin(), dl.end(), vec.begin());
    for( auto&& i : vec) {

        i = c;

    }

    vec[0] = make_cuDoubleComplex(0,0);

    //p(vec, dl);
    thrust::copy(vec.begin(), vec.end(), dl.begin());
    //f(du, vec);
    thrust::copy(du.begin(), du.end(), vec.begin());
    for( auto&& i : vec) {

        i = c;
    }

    vec[dl.size() - 1 ] = make_cuDoubleComplex(0, 0);

    //p(vec, du);
    thrust::copy(vec.begin(), vec.end(), du.begin());
}

void update_mdiag(thrust::device_vector<cuDoubleComplex>& d,
                  double tau,
                  double h,
                  double xmin) {


    cuDoubleComplex tauc = make_cuDoubleComplex( tau / 2.0, 0);
    cuDoubleComplex hconst = make_cuDoubleComplex(  1.0 / ( h * h), 0);

    thrust::host_vector<cuDoubleComplex> vec(d);

    //push_vec p();

    double x = xmin;

    cuDoubleComplex tmem = make_cuDoubleComplex(0,0);
    cuDoubleComplex one = make_cuDoubleComplex(1, 0);

    for( auto i = 0u; i < vec.size(); i++) {

        x += h * (double) i;
        tmem = tauc * (hconst + pot(x));
        tmem = make_cuDoubleComplex( -tmem.y, tmem.x);
        vec[i] = one + tmem;

    }
    //p(vec, d);
    thrust::copy(vec.begin(), vec.end(), d.begin());
} */
