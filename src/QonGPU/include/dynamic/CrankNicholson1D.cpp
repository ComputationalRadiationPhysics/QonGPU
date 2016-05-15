
#define DEBUG2(x) std::cout<<x<<std::endl

#include "CrankNicholson1D.hpp"
#include "CNKernels.h"




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
                                                  HDFile( 1),
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

void CrankNicholson1D::cusparse_init() {
    status = cusparseCreate( &handle);
    if(status != CUSPARSE_STATUS_SUCCESS) {
        std::cout<<"Error: Critical cuSparse unable to initialize"<<std::endl;
    }
    else {
        std::cout << "cuSparse successfully initialized!" << std::endl;
    }
}

void CrankNicholson1D::cusparse_destr() {
    // This function just destroys the cuSparse context
    // and throws an error if necessary

    status = cusparseDestroy( handle);
    if(status != CUSPARSE_STATUS_SUCCESS) {
        std::cout << "Error: cuSparse wasn't destroyed correctly!" << std::endl;
    }
    else {
        std::cout << "cuSparse was destroyed correctly!" << std::endl;
    }
}

void CrankNicholson1D::initfile(splash::DataCollector::FileCreationAttr& fa) {


    fa.fileAccType = splash::DataCollector::FAT_CREATE;
    HDFile.open(filename.c_str(), fa);
    splash::ColTypeDouble ctDouble;

    std::vector<double> p_sav(7);
    p_sav[0] = param->getxmax();
    p_sav[1] = param->getxmin();
    p_sav[2] = param->gettmax();
    p_sav[3] = param->getxmin();
    p_sav[4] = param->getnx();
    p_sav[5] = param->getnt();
    p_sav[6] = param->getz();
    splash::ColTypeDouble ctdouble;
    splash::Dimensions local(7, 1, 1);
    splash::Selection sel(local);
    HDFile.write(1,
                 ctdouble,
                 1,
                 sel,
                 "param_data",
                 p_sav.data());

}


void CrankNicholson1D::closefile() {
    HDFile.close();
}


void CrankNicholson1D::setstate(const thrust::host_vector<cuDoubleComplex>& v) {

    chunkl_d = v;
    chunkr_d = v;
    splash::DataCollector::FileCreationAttr fa;
    splash::DataCollector::initFileCreationAttr(fa);

    // Initialize the file and saves the paramerters
    initfile(fa);
    save_vectorh(0, v);
}

void CrankNicholson1D::savechunk(int step) {

    thrust::host_vector<cuDoubleComplex> vec_h(nx);
    // copy from Device to Host
    thrust::copy(chunkr_d.begin(), chunkr_d.end(), vec_h.begin());
    // Use STL Vectors to avoid raw pointer casts!
    std::vector<double> real(chunkr_d.size());
    std::vector<double> imag(chunkr_d.size());

    for(auto i = 0; i < chunkr_d.size(); ++i) {
        real[i] = vec_h[i].x;
        imag[i] = vec_h[i].y;
    }

    splash::ColTypeDouble type;
    splash::Dimensions dim(chunkr_d.size(),1,1);
    splash::Selection sel(dim);
    HDFile.write(step, type, 1, sel, "chunk_reals", real.data());
    HDFile.write(step, type, 1, sel, "chunk_img", imag.data());
}

void CrankNicholson1D::save_vector(int step,
                                   const thrust::device_vector<cuDoubleComplex>& v){
    thrust::host_vector<cuDoubleComplex> vec_h;
    // copy from Device to Host
    vec_h = v;
    std::vector<double> real(v.size());
    std::vector<double> imag(v.size());

    for(auto i = 0; i < v.size(); ++i) {
        real[i] = vec_h[i].x;
        imag[i] = vec_h[i].y;

    }

    splash::ColTypeDouble type;
    splash::Dimensions dim(v.size(),1,1);
    splash::Selection sel(dim);
    HDFile.write(step, type, 1, sel, "chunk_reals",real.data());
    HDFile.write(step, type, 1, sel, "chunk_img", imag.data());
}

void CrankNicholson1D::save_vectorh(int step, const thrust::host_vector<cuDoubleComplex>& v){

    std::vector<double> real(v.size());
    std::vector<double> imag(v.size());
    for(auto i = 0; i < v.size(); ++i) {
        real[i] = v[i].x;
        imag[i] = v[i].y;
    }

    splash::ColTypeDouble type;
    splash::Dimensions dim(v.size(),1,1);
    splash::Selection sel(dim);
    HDFile.write(step, type, 1, sel, "host_debug_reals",real.data());
    HDFile.write(step, type, 1, sel, "host_debug_img", imag.data());
}

void CrankNicholson1D::save_diag(int step, const thrust::device_vector<cuDoubleComplex>& v){
    thrust::host_vector<cuDoubleComplex> vec_h;
    // copy from Device to Host
    vec_h = v;
    std::vector<double> real(v.size());
    std::vector<double> imag(v.size());

    for(auto i = 0; i < v.size(); ++i) {
        real[i] = vec_h[i].x;
        imag[i] = vec_h[i].y;

    }

    splash::ColTypeDouble type;
    splash::Dimensions dim(v.size(),1,1);
    splash::Selection sel(dim);
    HDFile.write(step, type, 1, sel, "diag_reals",real.data());
    HDFile.write(step, type, 1, sel, "diag_img", imag.data());
}

void CrankNicholson1D::time_solve() {


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

    cusparse_init();
    for( auto i = 0; i < nt; ++i) {

        t += tau * (double) i;
        // check
        rhs_rt(c);
        // first perform the RHS Matrix multiplication!
        // Then update the non-constant main-diagonal!
        update_diagl<<< nx, 1>>>( dev_d, tau, h, xmin, nx);
        // right after that, we can call the cusparse Library
        // to write the Solution to the LHS chunk
        status2 = cusparseZgtsv( handle, nx, 1, dev_dl, dev_d, dev_du, dev_rhs, nx);

        DEBUG2("Currently calculation the "<<i<<"-th frame");
        thrust::copy( chunkr_d.begin(), chunkr_d.end(), chunkl_d.begin());

        savechunk(i+1);
    }
    cusparse_destr();
    closefile();
}
