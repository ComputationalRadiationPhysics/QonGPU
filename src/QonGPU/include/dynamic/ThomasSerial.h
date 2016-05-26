//
// Created by zaph0d on 25/05/16.
//

#pragma once

#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<vector>
#include<cassert>

#include"cuComplex.h"

using namespace thrust;
void solve_tridiagonal(const device_vector<cuDoubleComplex>& dev_du,
                       const device_vector<cuDoubleComplex>& dev_dl,
                       const device_vector<cuDoubleComplex>& dev_d,
                       device_vector<cuDoubleComplex>& dev_rhs) {

    host_vector<cuDoubleComplex> du(dev_du.size());
    host_vector<cuDoubleComplex> dl(dev_dl.size());
    host_vector<cuDoubleComplex> d(dev_d.size());
    host_vector<cuDoubleComplex> rhs(dev_rhs.size());
    host_vector<cuDoubleComplex> res(dev_rhs.size());


    thrust::copy( dev_du.begin(), dev_du.end(), du.begin());
    thrust::copy( dev_dl.begin(), dev_dl.end(), dl.begin());
    thrust::copy( dev_d.begin(), dev_d.end(), d.begin());
    thrust::copy( dev_rhs.begin(), dev_rhs.end(), rhs.begin());



    host_vector<cuDoubleComplex> du2(dev_du.size());
    host_vector<cuDoubleComplex> rhs2(dev_d.size());

    thrust::copy( du.begin(), du.end(), du2.begin());
    DEBUG2(rhs[100].x << " " << rhs[100].y);
    assert(abs(rhs[100].x) < 1e2);
    //thrust::copy( rhs.begin(), rhs.end(), rhs2.begin());

    du2[0] = du[0] / d[0];
    rhs2[0] = rhs[0] / d[0];

    cuDoubleComplex temp1 = make_cuDoubleComplex( 0, 0);
    cuDoubleComplex temp2 = make_cuDoubleComplex( 0, 0);

    // get maximum index
    const int m = du.size() -1 ;

    for(int i = 1; i < m; i++) {


        temp1 = d[i] - (dl[i] * du2[i - 1]);

        du2[i] = du[i] / temp1;

        temp1 = rhs[i] - (rhs2[ i -1 ] * dl[i]);
        //DEBUG2("temp1: "<< temp1.x << " " << temp1.y);
        temp2 = d[i] - (du2[i -1] * dl[i]);

        rhs2[i] = temp1 / temp2;

    }


    temp1 = rhs[m] - dl[m] * rhs2[m - 1];
    temp2 = d[m] - dl[m] * du2[m - 1];
    rhs2[m] = temp1 / temp2;

    res[m]  = rhs2[m];
    for(int i = m - 1; i >= 0; i--) {

        res[i] = rhs2[i] - du2[i] * res[i + 1];
    }

    thrust::copy(res.begin(), res.end(), dev_rhs.begin());
}