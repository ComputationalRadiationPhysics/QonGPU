// Currently under MIT license


#define TMAX 5000.0

#include "include/SimDef.hpp"




int main(int argc, char** argv) {

    if(argc < 2){
        std::cout << "Please specify filename!" << std::endl;
        std::cout << "Usage: qsolve <filename>" << std::endl;
        return -1;
    }
    std::string str = argv[1];

    double xma = 45.0;
    double xmi = -xma;
    double tma = TMAX;
    double tmi = 0.0;
    cudaDeviceReset();
    Params1D p(xma, xmi, tma, tmi, 1e5, 1e5, 1e8,1,str);
	SimDef<Numerov, CrankNicolson1D, Core1D, Params1D, 1> s(&p);
#ifndef TOTEST
    s.staticsolve();
    s.timerev();
#endif
    return EXIT_SUCCESS;

}








