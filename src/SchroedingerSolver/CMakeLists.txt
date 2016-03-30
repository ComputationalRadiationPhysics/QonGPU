######
#Right CMake Version
######

cmake_minimum_required(VERSION 3.0)

######
#Specify Project Name
######

project(qsolve)

######
#Set Compiler Flags
#####

set(CMAKE_CXX_FLAGS "--std=c++11")

######
#Find Packages
#####

#Find Cuda
find_package(CUDA REQUIRED)
set(LIBS ${LIBS} ${CUDA_LIBRARIES})
set(LIBS ${LIBS} ${CUDA_CUBLAS_LIBRARIES})
set(LIBS ${LIBS} ${CUDA_cusparse_LIBRARY})
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_20)

#Find MPI
#find_package(MPI REQUIRED)
#set(LIBS "${LIBS} ${MPI_LIBRARIES}")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS})

#Find HDF5 Package for file IO
FIND_PACKAGE(HDF5 COMPONENTS HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_HL_LIBRARIES})

#find Package Boost
find_package(Boost REQUIRED)
set(LIBS ${LIBS} ${BOOST_LIBRARIES})
include_directories(${BOOST_INCLUDE_DIR})

if(${CUDA_FOUND})
  message("-- Using CUDA")
  include_directories("${PROJECT_SOURCE_DIR}/include/")
  include_directories("${PROJECT_SOURCE_DIR}/include/*")
  link_directories("${PROJECT_SOURCE_DIR}/include/")
  file(GLOB SRC_FILES "main.cu")
  cuda_add_executable(qsolve ${SRC_FILES})
  target_link_libraries(qsolve ${LIBS} )

elseif(NOT ${CUDA_FOUND})

  message(FATAL_ERROR "CUDA NOT FOUND BUT IS REQUIRED")

 endif(${CUDA_FOUND})



