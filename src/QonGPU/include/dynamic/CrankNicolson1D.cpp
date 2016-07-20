
#define DEBUG2(x) std::cout<<x<<std::endl

//#define CUSPARSE_ON
//#define USE_SERIAL
#define USE_SPIKE
//#define MATRIX_OUTPUT
#include "CrankNicolson1D.hpp"


CrankNicolson1D::CrankNicolson1D(Params1D *_p):   param(_p),
                                                  nx( _p->getnx()),
                                                  nt( _p->getnt()),
                                                  E( 0.0),
                                                  chunkl_d( _p->getnx()),
                                                  chunkr_d( _p->getnx()),
                                                  tmax( _p->gettmax()),
                                                  tmin( _p->gettmin()),
                                                  xmax( _p->getxmax()),
                                                  xmin( _p->getxmin()),
                                                  d( _p->getnx()),
                                                  du( _p->getnx()),
                                                  dl( _p->getnx()),
                                                  filename( _p->getname()),
                                                  tau(( _p->gettmax() - _p->gettmin()) /_p->getnt())
{

}

CrankNicolson1D::~CrankNicolson1D() {}


void CrankNicolson1D::rhs_rt(double t) {

    // Prepare the rhs by using a triangular
    // matrix multiplication on rhs_rt
    // note that chunkl_d = chunkr_d since
    // only given chunkr_d to the routine would make
    // a temporal allocated field necessary!
	
    transform_rhs<<<512,3>>>(raw_pointer_cast(chunkl_d.data()),
                             raw_pointer_cast(chunkr_d.data()),
                             nx,
                             xmax,
                             xmin,
                             tau,
                             t);

}


void CrankNicolson1D::write_p(hid_t *f) {

    // Allocate and fill vector to save
    // Simulation parameters
    std::vector<double> p_sav(8);
    p_sav[0] = param->getxmax();
    p_sav[1] = param->getxmin();
    p_sav[2] = param->gettmax();
    p_sav[3] = param->gettmin();
    p_sav[4] = param->getnx();
    p_sav[5] = param->getnt();
    p_sav[6] = param->getz();
    p_sav[7] = param->geten();

    hsize_t rank = 1;
    hsize_t size = p_sav.size();
    H5LTmake_dataset(*f, "/params", rank, &size, H5T_NATIVE_DOUBLE,  p_sav.data() );

}


void saveblank(const thrust::device_vector<cuDoubleComplex>& v,
               hid_t *fl, int it) {

    thrust::host_vector<cuDoubleComplex> h(v.size());
    thrust::copy(v.begin(),v.end(), h.begin());
    std::vector<double> cs_rl(h.size());

    for(auto i = 0u; i < v.size(); ++i){

        cs_rl[i] = h[i].x;

    }

    std::string name = "/dset" + std::to_string(it) +"real";
    hsize_t rank = 1;
    hsize_t dim = cs_rl.size();
    H5LTmake_dataset(*fl, name.c_str(),rank, &dim,H5T_NATIVE_DOUBLE, cs_rl.data());

    for(auto i = 0u; i < v.size(); ++i){

        cs_rl[i] = h[i].y;

    }

    name = "/dset" + std::to_string(it) +"img";
    H5LTmake_dataset(*fl, name.c_str(),rank, &dim,H5T_NATIVE_DOUBLE, cs_rl.data());

}


void CrankNicolson1D::setstate(const thrust::host_vector<cuDoubleComplex>& v) {

    hid_t fl = H5Fcreate("copytest.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    thrust::copy(v.begin(), v.end(), chunkr_d.begin());
    saveblank(chunkr_d, &fl, 1);

    thrust::copy(v.begin(), v.end(), chunkl_d.begin());
    saveblank(chunkl_d, &fl, 0);
    H5Fclose(fl);
    /*
    thrust::host_vector<cuDoubleComplex> initial(nx);

    double h = (xmax - xmin)/nx;
    double ground = 0;
    double e = 0;
    DEBUG2("Setstate CALL");
    for(int i = 0; i < nx; i++) {

        e = xmin + i * h;
        ground = std::exp(-e*e/2)*std::pow(M_PI,0.25);
        chunkr_d[i] = make_cuDoubleComplex(ground, 0);
        chunkl_d[i] = make_cuDoubleComplex(ground, 0);
    }
    saveblank(chunkl_d, &fl, 0);
    H5Fclose(fl);*/
}


void CrankNicolson1D::time_solve() {

    // Define HDF5 File
    hid_t fl = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

#ifdef MATRIX_OUTPUT
    hid_t cfl = H5Fcreate("matr.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#endif

    // Write Parameters and define h as our discretization step
    write_p(&fl);
    double h = (xmax - xmin) / (double) nx;

    // constants of the diagonal
    const double c = - tau / (4.0 * h * h);

    // set t to starting point
    double t = param->gettmin();

    // cast to raw pointers for Kernel usage
    cuDoubleComplex * dev_d = raw_pointer_cast(d.data());
    cuDoubleComplex * dev_du = raw_pointer_cast(du.data());
    cuDoubleComplex * dev_dl = raw_pointer_cast(dl.data());
    cuDoubleComplex * dev_rhs = raw_pointer_cast(chunkr_d.data());


    // fill lower and upper diagonal
    // these stay constat for the whole time!

    create_const_diag <<< 512, 3 >>> (raw_pointer_cast(dl.data()),
                                      raw_pointer_cast(du.data()),
                                                                c,
                                                              nx);


	// create_diags(dl, du, make_cuDoubleComplex(0,c));
    // Save Diagonals for debug purposes
    //saveblank(dl, &cfl, 0);
    //saveblank(du, &cfl, 1);

    // Use cudaevent for time measurement!
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    float t_el = 0;

    // save copy!
    // saveblank(chunkl_d, &fl, 0);


#ifdef CUSPARSE_ON

    // Initialize cusparse if it's required
    cusparseStatus_t status;
    cusparseHandle_t handle = 0;
    status  = cusparseCreate(&handle);
    if(status != CUSPARSE_STATUS_SUCCESS) {

        std::cout<<"Error Init failed!"<<std::endl;

    }

#endif


    saveblank(chunkl_d, &fl, 0);
    for (int i = 0; i < nt; i++) {


        t = tau * (double) i + tmin;
		DEBUG2("Currently at t = "<< t);
        // Perform RHS multiplication
#ifdef MATRIX_OUTPUT
        saveblank(chunkr_d, &cfl, 2*i);
#endif
        rhs_rt(t-tau);
        //fast_mult(chunkr_d, tau, h, xmin);
#ifdef MATRIX_OUTPUT
        saveblank(chunkr_d, &cfl, 2*i+1);
#endif
        cuDoubleComplex check = chunkr_d[100];
        DEBUG2(check);



        // first perform the RHS Matrix multiplication!
        // Then update the non-constant main-diagonal!

        update_diagl <<< 512, 3 >>> (dev_d, tau, h, xmin, nx,  t);
        //update_mdiag(d, tau, h, xmin);
        //saveblank(d,  &cfl, i+1);
        //update_mdiag(d, tau, h, xmin);
        // right after that, we can call the cusparse Library
        // to write the Solution to the LHS chunkd
        cudaEventRecord(start);
        //fast_mult(chunkr_d, tau, h, xmin);
#ifdef USE_SPIKE

        gtsv_spike_partial_diag_pivot_v1<cuDoubleComplex, double>(dev_dl, dev_d, dev_du, dev_rhs, nx);
        DEBUG2("Spike Called!");

#endif

#ifdef CUSPARSE_ON

        status  = cusparseZgtsv(handle, nx, 1, dev_dl, dev_d, dev_du, dev_rhs, nx);
        if(status != CUSPARSE_STATUS_SUCCESS) {
            std::cout<<" Calulation went wrong"<<std::endl;
            assert(status == CUSPARSE_STATUS_SUCCESS);
        }

#endif


#ifdef USE_SERIAL

        solve_tridiagonal(du, dl, d, chunkr_d);

#endif

        cudaEventRecord(stop);

        // Write RHS to LHS
        thrust::copy(chunkr_d.begin(), chunkr_d.end(), chunkl_d.begin());


        // Calculate Time
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&t_el, start, stop);

        // Debug Messages
        std::cout << "Generated the " << i << "-th frame" << std::endl;
        std::cout << "Frame generation time: " << t_el << "ms" << std::endl;


        assert(check.x < 100);
        assert(check.y < 100);

        if (i % 1000 == 0)
            saveblank(chunkr_d, &fl, i + 1);


        if (i == 1e9) {

            saveblank(chunkr_d, &fl, 1e5);
            i = 2 * nt;

        }

        //saveblank(chunkl_d, &fl, i + 1);
    }

    std::cout << "The starting Energy was: " << param->geten() << std::endl;

#ifdef CUSPARSE_ON
    // Destroy cusparse
    status = cusparseDestroy(handle);
    if(status != CUSPARSE_STATUS_SUCCESS){

        std::cout<<"Error cusparse couldn't be destroyed!"<<std::endl;

    }
#endif

    H5Fclose(fl);

#ifdef MATRIX_OUTPUT

    H5Fclose(cfl);

#endif
    DEBUG2("h was: "<<h);
    DEBUG2("tau was: "<< tau);
    DEBUG2("C was: "<< c);
    DEBUG2("Z : "<< param->getz());

}
