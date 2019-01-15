cmake_minimum_required(VERSION 3.10)
project(QonGPU)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/QonGPU/include)
include_directories(src/QonGPU/include/dynamic)
include_directories(src/QonGPU/include/output)
include_directories(src/QonGPU/include/params)
include_directories(src/QonGPU/include/potentials)
include_directories(src/QonGPU/include/static)

add_executable(QonGPU
        src/QonGPU/include/dynamic/CNKernels.h
        src/QonGPU/include/dynamic/ComplexOperators.h
        src/QonGPU/include/dynamic/CrankNicolson1D.cpp
        src/QonGPU/include/dynamic/CrankNicolson1D.hpp
        src/QonGPU/include/dynamic/cusparse_ops.hxx
        src/QonGPU/include/dynamic/MatrixGeneration.h
        src/QonGPU/include/dynamic/MemFunctors.h
        src/QonGPU/include/dynamic/spike_kernel.hxx
        src/QonGPU/include/dynamic/ThomasSerial.h
        src/QonGPU/include/dynamic/TimeOperator.hpp
        src/QonGPU/include/dynamic/TimeOperator1D.hpp
        src/QonGPU/include/dynamic/TimeOperator2D.hpp
        src/QonGPU/include/dynamic/TimeOperator3D.hpp
        src/QonGPU/include/dynamic/TridiagMult.h
        src/QonGPU/include/output/IOHandle.hpp
        src/QonGPU/include/output/IOHandle1D.cpp
        src/QonGPU/include/output/IOHandle1D.h
        src/QonGPU/include/params/Params.hpp
        src/QonGPU/include/params/Params1D.hpp
        src/QonGPU/include/params/Params2D.hpp
        src/QonGPU/include/params/Params3D.hpp
        src/QonGPU/include/potentials/Core1D.hpp
        src/QonGPU/include/potentials/Potential.hpp
        src/QonGPU/include/potentials/Potential1D.hpp
        src/QonGPU/include/potentials/Potential2D.hpp
        src/QonGPU/include/potentials/Potential3D.hpp
        src/QonGPU/include/static/Numerov.cpp
        src/QonGPU/include/static/Numerov.hpp
        src/QonGPU/include/static/Numerov1D.hpp
        src/QonGPU/include/static/StaticSolver.hpp
        src/QonGPU/include/static/StaticSolver1D.hpp
        src/QonGPU/include/static/StaticSolver2D.hpp
        src/QonGPU/include/static/StaticSolver3D.hpp
        src/QonGPU/include/AllHeader.hpp
        src/QonGPU/include/Domain.hpp
        src/QonGPU/include/SimDef.hpp
        CMakeLists.txt)
