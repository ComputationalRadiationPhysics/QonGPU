#include <vector>
#include <iostream>
#include "Numerov.hpp"



Numerov::Numerov():nx(0),ne(0),xmax(0),xmin(0){}

Numerov::Numerov(Params1D *pa,complex<double>* ps): param(pa),
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
    // Create the Device Pointer calculating the chunks
    double* dev_ptr;
    int dev_ne = 0;
    // Next let's allocate the required chunk memory on the device side
    cudaMalloc((void**)&dev_ptr,sizeof(double)*nx*CHUNKSIZE);
    // Make use of some local variables
    int index = 0;
    double dE = V(0,0,z)/ (double)ne;
    // This will be the starting energy for each chunk calculation
    double En = 0;
    while( index < ne) {
        //copy initals on device
        dev_ne = CHUNKSIZE;
        if( index + CHUNKSIZE > ne){
            dev_ne = ne - index;
            cudaFree(dev_ptr);
            cudaMalloc((void**) &dev_ptr, sizeof(double) * nx * dev_ne);
        }
        DEBUG("Calculating chunk: "<< index/CHUNKSIZE)
        cudaMemcpy( dev_ptr, cache.data(), sizeof(double)*nx*dev_ne, cudaMemcpyHostToDevice);
        En = dE * (double) index;
        DEBUG("Calculating with starting energy: " << En)
        iter1<<< dev_ne, 1>>>( dev_ptr, nx, dev_ne, xmax, xmin, z, En, dE);
        cudaMemcpy( chunk.data(), dev_ptr, sizeof(double)*nx*dev_ne, cudaMemcpyDeviceToHost);
        bisect(index);
        index += CHUNKSIZE;
    }
    savelevels();
}

void Numerov::savelevels(){
    // Function to provide saving functionality of the energy Levels
    // First Allocate an appropiate vector
    vector<double> buffer1(results.size()*nx);
    vector<double> buffer2(eval.size());
    vector<double> buffer3(nx);
    // Save the results into buffer1 vector
    DEBUG(results.size())
    results.pop_front();
    for(auto i = 0; i < results.size();i+=nx){
        // Write the first list elements into a vector
        for(auto i = 0; i < nx; i++) {
            buffer3[i] = results.front()[i];
        }
        // Now we write the data into another vector
        // I know this is stupid, but since list has no direct data access
        // I see no other choice.
        for(auto j = 0; j < nx; j++) {
            buffer1[i+j] = buffer3[j];
        }
        // Now pop the first element from the list
        results.pop_back();
    }
    for(auto it = 0; it < eval.size(); it++) {
        // We do the analog thing for the enegy
        // It's a lot simpler!
        buffer2[it] = eval.back();
        eval.pop_back();
    }
    // Create a new HDF5 file
    hid_t file_id;
    file_id = H5Fcreate("static_results.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    hsize_t dims = buffer1.size();
    // Create a HDF5 Data set and write buffer1
    H5LTmake_dataset(file_id, "/numres", 1, &dims, H5T_NATIVE_DOUBLE, buffer1.data());
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
}

bool Numerov::sign(double s){
    return std::signbit(s);
}



void Numerov::bisect(int j) {
    // Iterate through chunk data
    // create local variable for the offset
    // off is the index of the last Element
    DEBUG("Bisecting")
    const int off = nx - 1;
    auto it = chunk.begin();
    vector<double> temp( nx);
    for( auto i = 2 * off; i < chunk.size(); i += nx){
        if(sign( chunk[ i ]) != sign( chunk[ i - nx])){
            DEBUG("Signchange deteceted!")
            if( (fabs( chunk[ i]) < fabs( chunk[i - nx]))&&chunk[i]<1e-3) {
                std::copy( it + i, it + i + nx, temp.begin());
                results.push_back( temp);
                // @TODO energy output
                std::cout << "Energy level found" << std::endl;
            }
            else {
                if(chunk[i-nx]<1e-3) {
                    std::copy(it + i, it + i + nx, temp.begin());
                    results.push_back(temp);
                    std::cout << "Energy level found" << std::endl;
                }

            }
        }
    }
}
void Numerov::tempprint(){
    hid_t fileid;
    hsize_t dim=7;
    //Create temporal array for parameter location:
    double* tempdata = (double*) malloc(sizeof(double)*7);
    tempdata[0]=param->getxmax();
    tempdata[1]=param->getxmin();
    tempdata[2]=param->gettmax();
    tempdata[3]=param->gettmin();
    tempdata[4]=(double) param->getnx();
    tempdata[5]=(double) param->getnt();
    tempdata[6]=(double) param->getne();
    //create the hdf5 file:
    fileid = H5Fcreate("res_temp.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    //Print the parameters in the file
    H5LTmake_dataset(fileid,"/params1d",1,&dim, H5T_NATIVE_DOUBLE,tempdata);
    free(tempdata);
    //rewrite dims
    dim = (param->getnx())*(CHUNKSIZE);
    //print the simulation results in the HDF5 file
    H5LTmake_dataset(fileid, "/numres", 1, &dim, H5T_NATIVE_DOUBLE,chunk.data());
    //now write the energy indices to the file
    //close the file
    H5Fclose(fileid);
    //free memory
    DEBUG("Finished saving!")
};