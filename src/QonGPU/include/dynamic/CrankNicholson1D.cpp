
#define DEBUG2(x) std::cout<<x<<std::endl

#include "CrankNicholson1D.hpp"
#include "CNKernels.h"
#include "hdf5.h"
#include "hdf5_hl.h"



CrankNicholson1D::CrankNicholson1D(Params1D *_p): param(_p),
                                                  nx( _p->getnx()),
                                                  nt( _p->getnt()),
                                                  E( 0.0),
                                                  chunk_h(_p->getnx()),
                                                  chunkr_d( _p->getnx()),
                                                  tmax( _p->gettmax()),
                                                  tmin( _p->gettmin()),
                                                  xmax( _p->getxmax()),
                                                  xmin( _p->getxmin()),
                                                  d( _p->getnx()),
                                                  du( _p->getnx()),
                                                  dl( _p->getnx()),
                                                  inital( _p->getnx()),
                                                  filename( _p->getname()),
                                                  tau(( _p->gettmax()- _p->gettmin()) /_p->getnt())
{

}

CrankNicholson1D::~CrankNicholson1D() {}

void CrankNicholson1D::rhs_rt( const double c) {
    // Prepare the rhs by using a triangular
    // matrix multiplication on rhs_rt
    // note that chunkl_d = chunkr_d since
    // only given chunkr_d to the routine would make
    // a temporal allocated field necessary!

    transform_rhs<<<nx,1>>>(raw_pointer_cast(chunkl_d.data()),
            raw_pointer_cast(chunkr_d.data()),
            nx, xmax,
            xmin,tau,c);
}


void CrankNicholson1D::write_p(hid_t *f) {




    std::vector<double> p_sav(7);
    p_sav[0] = param->getxmax();
    p_sav[1] = param->getxmin();
    p_sav[2] = param->gettmax();
    p_sav[3] = param->getxmin();
    p_sav[4] = param->getnx();
    p_sav[5] = param->getnt();
    p_sav[6] = param->getz();

}

void CrankNicholson1D::setstate(const thrust::host_vector<cuDoubleComplex>& v) {

    chunkl_d = v;
    chunkr_d = v;
}

void saveblank(const thrust::device_vector<cuDoubleComplex>& v,
               hid_t *fl, int it){

    thrust::host_vector<cuDoubleComplex> h = v;
    std::vector<double> cs_rl(h.size());
    for(auto i = 0; i < v.size(); ++i){
        cs_rl[i] = h[i].x;
    }
    std::string name = "/dset" + std::to_string(it) +"real";
    hsize_t rank = 1;
    hsize_t dim = cs_rl.size();
    H5LTmake_dataset(*fl, name.c_str(),rank, &dim,H5T_NATIVE_DOUBLE, cs_rl.data());
    for(auto i = 0; i < v.size(); ++i){
        cs_rl[i] = h[i].y;
    }
    name = "/dset" + std::to_string(it) +"img";
    H5LTmake_dataset(*fl, name.c_str(),rank, &dim,H5T_NATIVE_DOUBLE, cs_rl.data());

}
void CrankNicholson1D::time_solve() {


    hid_t fl = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);



    // This routine is now slightly longer

    const double h = (xmax - xmin) / (double) nx;
    // constants of the diagonal
    const double c =  - tau / (4.0 * pow(h,2.0));
    DEBUG2("h equals " << h);
    DEBUG2("pow(h,2) = "<< pow(h,2.0));
    DEBUG2("Tau Equals " << tau);
    DEBUG2("C = "<< c);
    double t = param->gettmin();

    //splash::DataCollector::FileCreationAttr fa;
    //splash::DataCollector::initFileCreationAttr(fa);

    // Initialize the file and saves the paramerters
    //initfile(fa);

    cuDoubleComplex* dev_d = raw_pointer_cast(d.data());
    cuDoubleComplex* dev_du = raw_pointer_cast(du.data());
    cuDoubleComplex* dev_dl = raw_pointer_cast(dl.data());
    cuDoubleComplex* dev_rhs = raw_pointer_cast(chunkr_d.data());


    // Check
    create_const_diag<<<nx ,1>>>( raw_pointer_cast(dl.data()),
            raw_pointer_cast(du.data()),
            c , nx);
    //saveblank(chunkl_d, &fl, 0);
    for( auto i = 0; i < nt; ++i) {

        t += tau * (double) i;
        // check
        rhs_rt(c);
        // first perform the RHS Matrix multiplication!
        // Then update the non-constant main-diagonal!
        update_diagl<<< nx, 1>>>( dev_d, tau, h, xmin, nx);
        // right after that, we can call the cusparse Library
        // to write the Solution to the LHS chunkd
        gtsv_spike_partial_diag_pivot_v1<cuDoubleComplex,double>(dev_dl, dev_d, dev_du, dev_rhs,nx);

        chunkl_d = chunkr_d;
        saveblank(chunkl_d, &fl, i);

    }
    H5Fclose(fl);
}
