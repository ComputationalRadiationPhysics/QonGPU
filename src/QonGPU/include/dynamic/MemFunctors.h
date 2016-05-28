/*
 * These Function Objects are an generalized interface
 * for the CUDA host/Device memory model. We wan't to use both
 * thrust vectors as well as simple CUDA Device pointers, therefore
 * we pack everything in a fetch and push functor!
 *
 */
#pragma once

template<typename DC, typename HC,typename T>
struct fetch {

    void operator()(const DC& in, HC& out) {
        std::cout<<"Warning! No specific operation specified, now killing process!"<<std::endl;
        assert(1==2);
    };
};


template<>
struct fetch<thrust::device_vector<cuDoubleComplex> ,
             thrust::host_vector<cuDoubleComplex>,
             cuDoubleComplex >
{


    void operator()(const thrust::device_vector<cuDoubleComplex>& in,
                    thrust::host_vector<cuDoubleComplex>& out) {

        thrust::copy(in.begin(), in.end(), out.begin());

    }

};

typedef  fetch<thrust::device_vector<cuDoubleComplex> ,
         thrust::host_vector<cuDoubleComplex>,
         cuDoubleComplex > fetch_vec;

template<typename DC, typename HC, typename T>
struct push {

    void operator()(const HC& in, DC& out) {

        std::cout<<"Warning! No specific operation specified, now killing process!"<<std::endl;
        assert(1==2);

    }
};


template<>
struct push<thrust::device_vector<cuDoubleComplex>,
           thrust::host_vector<cuDoubleComplex>,
            cuDoubleComplex>
{

        void operator()(const thrust::host_vector<cuDoubleComplex>& in,
                        thrust::device_vector<cuDoubleComplex>& out) {

            thrust::copy(in.begin(), in.end(), out.begin());
        }
};

typedef push<thrust::device_vector<cuDoubleComplex>,
             thrust::host_vector<cuDoubleComplex>,
             cuDoubleComplex> push_vec;


