//
// Created by zaph0d on 28/05/16.
//
#pragma  once

void fast_mult( thrust::device_vector<cuDoubleComplex> psi_d,
                double tau, double h, double xmin) {

    fetch_vec f;
    push_vec p;

    thrust::host_vector<cuDoubleComplex> psi_h(psi_d);
    thrust::host_vector<cuDoubleComplex> cpy(psi_d.size(), make_cuDoubleComplex(0,0));

    int n = psi_h.size();
    double x = xmin;
    double xmax = xmin + h * psi_h.size();
    const cuDoubleComplex hconst = make_cuDoubleComplex( -1.0 / ( 2.0 * h * h), 0);
    const cuDoubleComplex tauc = make_cuDoubleComplex( -tau / 2.0, 0);

    cuDoubleComplex helper = make_cuDoubleComplex(0, 0);
    cuDoubleComplex helper2 = make_cuDoubleComplex(0, 0);

    helper = hconst * tauc * ( psi_h[1] - 2 * psi_h[0]);
    helper2 = tauc * pot(0) * psi_h[0];
    helper = helper + helper2;
    helper = make_cuDoubleComplex(-helper.y, helper.x);
    cpy[0] = psi_h[0] + helper;

    for(int i = 1; i < psi_h.size(); i++) {

        x += h * i;
        helper = hconst * tauc * ( psi_h[i+1] + psi_h[i-1] - psi_h[i]);
        helper2 = tauc * pot(x) * psi_h[i];
        helper = helper + helper2;
        helper = make_cuDoubleComplex( -helper.y, helper.x);
        cpy[i] = psi_h[i] + helper;

    }

    helper = hconst * tauc * ( psi_h[n-2] - 2*psi_h[n-1]);
    helper2 = tauc * pot(xmax) * psi_h[n-1];
    helper = helper + helper2;
    helper = make_cuDoubleComplex( -helper.y, helper.x);
    cpy[n-1] = psi_h[n - 1] + helper;

    p(cpy, psi_d);



}