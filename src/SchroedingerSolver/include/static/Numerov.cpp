#include <vector>
#include <iostream>
#include "Numerov.hpp"


Numerov::Numerov():nx(0),ne(0),xmax(0),xmin(0){}

Numerov::Numerov(Params1D *pa,std::complex<double>* ps): param(pa),
                                                    chunk(pa->getnx()*CHUNKSIZE),cache(pa->getnx()*CHUNKSIZE),
                                                    nx(pa->getnx()),
                                                    ne(pa->getne()),
                                                    z(pa->getz()),
                                                    xmax(pa->getxmax()),
                                                    xmin(pa->getxmin()) {
    // initialize the cache, with the inital values
    for(auto it = cache.begin(); it != cache.end(); it += nx) {
        // the first entry always has to be zero, the second
        // is an arbitrary small value
        *(it) = 0;
        *(it+1) = -1e-10;
    }
    //get the minimal energy E_n is in [V(0,0,z),0]
    //We'll define it as positive
    //Kernel code then uses negaive values!
    E = V(0,0,z);
}

Numerov::~Numerov(){}

void Numerov::solve(){
    // This Loop is used to create
    for(auto j = 1; j < 252; j++) {
        z = j;
        DEBUG("Solving for Z ="<<z)

        double *dev_ptr;

        int dev_ne = 0;
        // Next let's allocate the required chunk memory on the device side
        cudaMalloc((void **) &dev_ptr, sizeof(double) * nx * CHUNKSIZE);
        // Make use of some local variables
        int index = 0;
        double dE = V(0, 0, z) / (double) ne;
        // This will be the starting energy for each chunk calculation
        double En = 0;
        while (index < ne) {
            //copy initals on device
            dev_ne = CHUNKSIZE;
            if (index + CHUNKSIZE > ne) {
                dev_ne = ne - index;
                cudaFree(dev_ptr);
                cudaMalloc((void **) &dev_ptr, sizeof(double) * nx * dev_ne);
            }
            DEBUG("Calculating chunk: " << index / CHUNKSIZE)
            cudaMemcpy(dev_ptr, cache.data(), sizeof(double) * nx * dev_ne, cudaMemcpyHostToDevice);
            En = dE * (double) index;
            DEBUG("Calculating with starting energy: " << En)
            iter1 << < dev_ne, 1 >> > (dev_ptr, nx, dev_ne, xmax, xmin, z, En, dE);
            cudaMemcpy(chunk.data(), dev_ptr, sizeof(double) * nx * dev_ne, cudaMemcpyDeviceToHost);
            bisect(En);
            index += CHUNKSIZE;
        }
    }
    // After all the calculations done we can save our energy levels!
    savelevels();
}

void Numerov::savelevels(){
    // Function to provide saving functionality of the energy Levels
    // First Allocate an appropiate vector

    DEBUG(res.size()/nx)

    hid_t file_id;
    vector<double> buffer2(eval.size());


    for(auto it = 0; it < eval.size(); it++) {
        // We do the analog thing for the enegy
        // It's a lot simpler!
        buffer2[it] = eval.back();
        eval.pop_back();
    }
    // Create a new HDF5 file

    file_id = H5Fcreate("static_results.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    hsize_t dims = res.size();

    // Create a HDF5 Data set and write buffer1
    H5LTmake_dataset(file_id, "/numres", 1, &dims, H5T_NATIVE_DOUBLE, res.data());
    // Analog for buffer2
    dims = buffer2.size();
    H5LTmake_dataset(file_id, "/evals", 1, &dims, H5T_NATIVE_DOUBLE, buffer2.data());
    //Save some necessary parameters
    vector<double> p(3);
    p[0] = nx;
    p[1] = ne;
    p[2] = xmax;
    dims = 3 ;
    H5LTmake_dataset( file_id, "/params", 1, &dims, H5T_NATIVE_DOUBLE, p.data());
    // close the file
    H5Fclose(file_id);
    DEBUG("CALL END")
}


bool Numerov::sign(double s){
    return std::signbit(s);
}



void Numerov::bisect(double j) {
    // Iterate through chunk data
    // create local variable for the offset
    // off is the index of the last Element
    const int off = nx - 1;
    auto it = chunk.begin();
    vector<double> temp( nx);
    double dE = -V(0,0,z)/ne;
    double En = 0;
    for( auto i = 2 * off; i < chunk.size(); i += nx){
        if(sign( chunk[ i ]) != sign( chunk[ i - nx])){
            if( (fabs( chunk[ i]) < fabs( chunk[i - nx]))&&chunk[i]<1e-3) {
                std::copy( it + i, it + i + nx, temp.begin());
                results.push_back( temp);
                std::cout << "Energy level found" << std::endl;
                std::cout << "Detected energy level: "<< En << std::endl;
                res.resize(res.size()+nx);
                auto iter = res.end()-nx;
                std::copy(it+i,it+i+nx,iter);
                En = ( -j + i/nx * dE);
                eval.push_back(En);
            }
            else {
                if(chunk[i-nx]<1e-3) {
                    std::copy(it + i, it + i + nx, temp.begin());
                    results.push_back(temp);
                    results.pop_back();
                    res.resize(res.size()+nx);
                    auto iter = res.end()-nx;
                    std::copy(it+i,it+i+nx,iter);
                    std::cout << "Energy level found" << std::endl;
                    std::cout << "Detected energy level: "<< En << std::endl;
                    En = (-j + (i/nx - 1) * dE);
                    eval.push_back(En);
                }
            }
        }
    }
}
