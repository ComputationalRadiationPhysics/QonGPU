
#ifndef ALLHEADER_H
#define ALLHEADER_H


#include <iostream>
#include <complex>

#include"params/Params.hpp"
#include "params/Params1D.hpp"
#include "params/Params2D.hpp"
#include "params/Params3D.hpp"

#include "dynamic/TimeOperator.hpp"
#include "dynamic/TimeOperator1D.hpp"
#include "dynamic/TimeOperator2D.hpp"
#include "dynamic/TimeOperator3D.hpp"
#include "dynamic/CrankNicolson1D.cpp"
#include "dynamic/ComplexOperators.h"

#include "static/StaticSolver.hpp"
#include "static/StaticSolver1D.hpp"
#include "static/StaticSolver2D.hpp"
#include "static/StaticSolver3D.hpp"

#include "potentials/Potential.hpp"
#include "potentials/Potential1D.hpp"
#include "potentials/Potential2D.hpp"
#include "potentials/Potential3D.hpp"
#include "potentials/Core1D.hpp"

#include "output/IOHandle.hpp"
#include "output/IOHandle1D.h"
#include "output/IOHandle1D.cpp"
//#include "static/Numerov1D.hpp"
#include "static/Numerov.cpp"

#endif
