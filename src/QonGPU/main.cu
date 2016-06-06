// Currently under MIT license

//#define DEBUG(x) std::cout<<x<<std::endl;
//#define STATUS(x) std::cout<<x<<"...";
//#define ENDSTATUS std::cout<<"DONE!"<<std::endl;

#include "include/SimDef.hpp"


int main(int argc, char** argv) {

    if(argc < 2){
        std::cout << "Please specify filename!" << std::endl;
        std::cout << "Usage: qsolve <filename>" << std::endl;
        return -1;
    }
    std::string str = argv[1];

    double xma = 4000.0;
    double xmi = -4000.0;
    double tma = 3190.0;
    double tmi = 0.0;
    cudaDeviceReset();
    Params1D p(xma, xmi, tma, tmi,(size_t) 3.2*1e6, (size_t) 1.45*1e6, 1e8,1,str);

    SimDef<Numerov, CrankNicolson1D, Core1D, Params1D, 1> s(&p);
#ifndef TOTEST
    DEBUG2("CALL 2");
    s.staticsolve();
    DEBUG2("CALL 3");
    s.timerev();
#endif
    return EXIT_SUCCESS;

}








