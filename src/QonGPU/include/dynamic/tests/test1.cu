//
// Created by max on 29/04/16.
//
#define DEBUG(x) std::cout<<x<<std::endl;
#define STATUS(x) std::cout<<x<<"...";
#define ENDSTATUS std::cout<<"DONE!"<<std::endl;

#include <iostream>
#include <complex>
#include <vector>



#include "../../params/Params1D.hpp"
#include "../TimeOperator1D.hpp"
#include "../TimeOperator.hpp"
#include "../CrankNicolson1D.cpp"

#include "../CNKernels.h"

#define BOOST_TEST_MODULE "DynTest"
#include "boost/test/included/unit_test.hpp"



/*
BOOST_AUTO_TEST_CASE(devicefunctions) {
// check the potential function
    cuDoubleComplex c0 = make_cuDoubleComplex(1.0,0);
    BOOST_CHECK_EQUAL(pot(0).x,c0.x);
    BOOST_CHECK_EQUAL(pot(0).y,c0.y);
// Check the outputs of the mult_rhs method
    const double con1 = 2.0;
    const double con2 = 3.0;
    const double con3 = 4.0;
    const double con4 = 0.0;
    const double con5 = 6.0;
    const double con6 = 0.0;
    cuDoubleComplex c1 = make_cuDoubleComplex(2.0, 2.0);
    cuDoubleComplex c2 = make_cuDoubleComplex(3.0, 3.0);
    cuDoubleComplex c3 = make_cuDoubleComplex(4.0, 4.0);
    cuDoubleComplex c4 = make_cuDoubleComplex(con4, con4);
    cuDoubleComplex c5 = make_cuDoubleComplex(con5, con5);
    cuDoubleComplex c6 = make_cuDoubleComplex(con6, con5);
    cuDoubleComplex h1 = make_cuDoubleComplex(1.0 ,0);
    cuDoubleComplex h2 = make_cuDoubleComplex(1.0, 0);
    mult_rhs( &c3, &c2, &c1, &c5, &c4, &c6, h1, h2 ,0);
    BOOST_CHECK_EQUAL(c5.x,3.0);
    BOOST_CHECK_EQUAL(c5.y,5.0);

    // checking the transform RHS
    cuDoubleComplex s1  = make_cuDoubleComplex(1.0 ,0);
    cuDoubleComplex s2 = make_cuDoubleComplex(3.0,0);
    cuDoubleComplex d;
    const double c_rhs = 5.0;
    const double x_rhs = 0;
    cuDoubleComplex t1 = make_cuDoubleComplex( 1.0, 2.0);
    transform_diag( d, s1, s2, c_rhs, x_rhs,t1);
    BOOST_CHECK_EQUAL(d.x,19.0);
    BOOST_CHECK_EQUAL(d.y,-9.0);
}
*/

BOOST_AUTO_TEST_CASE(BinaryDeviceOperators) {
        cuDoubleComplex res = make_cuDoubleComplex(4.5,1.0);
        cuDoubleComplex l1 = make_cuDoubleComplex(1.0,1.0);
        cuDoubleComplex l2 = make_cuDoubleComplex(2.0,1.0);

        res = l1 + l2 + res;

        BOOST_CHECK_EQUAL(res.x, 7.5);
        BOOST_CHECK_EQUAL(res.y, 3.0);

        std::cout<< res.x<<" "<<res.y<<std::endl;

        res = make_cuDoubleComplex(2.0,1.9);
        l1 = make_cuDoubleComplex(4.0,5.0);
        l2 = make_cuDoubleComplex(5.0,5.0);


        res = res*l1*l2;
        BOOST_CHECK_EQUAL(res.x, -95.5);
        BOOST_CHECK_EQUAL(res.y, 80.5);
        std::cout<< res.x<<" "<<res.y<<std::endl;

        res = make_cuDoubleComplex(2.0,1.9);
        l1 = make_cuDoubleComplex(4.0,5.0);
        res = res / l1;
        BOOST_REQUIRE_CLOSE(res.x, 0.426829268292683, 1e-4);
        BOOST_REQUIRE_CLOSE(res.y, -0.05853658536585366, 1e-4);

        std::cout<< res.x<<" "<<res.y<<std::endl;

        res = make_cuDoubleComplex(2.0,1.9);
        l1 = make_cuDoubleComplex(4.0,5.0);

        res = res - l1;
        BOOST_REQUIRE_CLOSE(res.x, -2.0, 1e-4);
        BOOST_REQUIRE_CLOSE(res.y, -3.1, 1e-4);
}






