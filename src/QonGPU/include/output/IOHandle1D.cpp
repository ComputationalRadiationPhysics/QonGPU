//
// Created by zaph0d on 21/04/16.
//

#include "IOHandle1D.h"



IOHandle1D::IOHandle1D(Params1D *_p) : real_cache( _p->getnx()),
                                       img_cache( _p->getnx()){
    /*
     * Create the HDF5 file, during object life-time!
     */

}