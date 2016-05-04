#pragma  once

#include <thrust/host_vector.h>
#include <cuComplex.h>

#include "TimeOperator.hpp"



class TimeOperator1D : TimeOperator {
public:
    virtual void setstate(const thrust::host_vector<cuDoubleComplex>& v) = 0;
};
