//
// Created by zaph0d on 21/04/16.
//

#include "IOHandle1D.h"


//First let's define the constructors
IOHandle1D::IOHandle1D(): nx(0),ne(0),nt(0){ }
IOHandle1D::IOHandle1D(Params1D *p):nx(p->getnx()),
                                    ne(p->getne()),
                                    nt(p->getnt()),
                                    param(p),
                                    chunk(p->getnx())
{
    initfile();
}
IOHandle1D::~IOHandle1D() {
    closefile();
}
//Open and close wrappers for file!
//Default lifetime of file is the lifetime of IOHandle1D
void IOHandle1D::initfile() {
    file = H5Fcreate(FILENAME,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
}
void IOHandle1D::closefile() {
    H5Fclose(file);
}

void IOHandle1D::copychunk(const std::vector<std::complex<double>>& c){
    // Create a iterator for loop to get a efficient as possible chunk
    // copy
    for(auto i =  c.begin(); i != c.end(); ++i) {
        auto d = std::distance(c.begin(),i);
        chunk.r[d] = i->real();
        chunk.i[d] = i->imag();
    }
}
void IOHandle1D::savechunk(const std::vector<std::complex<double>>& c) {
    // First let's do a copy operation to an struct of arrays
    copychunk(c);

}

