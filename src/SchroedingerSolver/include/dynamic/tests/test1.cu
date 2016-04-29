//
// Created by max on 29/04/16.
//
#define DEBUG(x) std::cout<<x<<std::endl;
#define STATUS(x) std::cout<<x<<"...";
#define ENDSTATUS std::cout<<"DONE!"<<std::endl;

#include <iostream>
#include <complex>
#include <vector>



#include <params/Params1D.hpp>
#include "../TimeOperator1D.hpp"
#include "../TimeOperator.hpp"
#include "../CrankNicholson1D.cpp"

#include <dynamic/CNKernels.h>

#define BOOST_TEST_MODULE "DynTest"
#include "boost/test/included/unit_test.hpp"


BOOST_AUTO_TEST_CASE(constructor) {

        std::complex<double> xma = 0.0;
    std::complex<double> xmi = 0.0;
    std::complex<double> tmi = 0.0;
    std::complex<double> tma = 0.0;
    const  size_t s1 = 3;
    const size_t s2 = 5;
    const size_t s3 = 6;
    const size_t s4 = 2;
    vector<cuDoubleComplex> v(s1);
    Params1D p( xma, xmi, tma, tmi, s1, s2, s3, s4);
    CrankNicholson1D cn(&p, v);
    BOOST_CHECK_EQUAL(cn.getnx(), s1);
    BOOST_CHECK_EQUAL(cn.getnt(), s3);
    BOOST_CHECK_EQUAL(cn.getxmax(),0.0);
    BOOST_CHECK_EQUAL(cn.getxmin(),0.0);
    BOOST_CHECK_EQUAL(cn.gettmax(),0.0);
    BOOST_CHECK_EQUAL(cn.gettmin(),0.0);
    cn.cusparse_init();
    cn.cusparse_destr();


}

BOOST_AUTO_TEST_CASE(devicefunctions) {
    cuDoubleComplex c0 = make_cuDoubleComplex(1.0,0);
    BOOST_CHECK_EQUAL(pot(0).x,c1.x);
    BOOST_CHECK_EQUAL(pot(0).y,c1.y);

}


void dummy(){
    const double con1 = 2.0;
    const double con2 = 3.0;
    const double con3 = 4.0;
    const double con4 = 5.0;
    const double con5 = 6.0;
    const double con6 = 7.0;
    cuDoubleComplex c1 = make_cuDoubleComplex(con1, con1);
    cuDoubleComplex c2 = make_cuDoubleComplex(con2, con2);
    cuDoubleComplex c3 = make_cuDoubleComplex(con3, con3);
    cuDoubleComplex c4 = make_cuDoubleComplex(con4, con4);

}






