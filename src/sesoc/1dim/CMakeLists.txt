#######################
#Required CMake Version
#######################
cmake_minimum_required(VERSION 3.0)

###################
#Add Subdirectories
###################


##################
#Project definitions
##################

project(1dim)

##################
#compiler flag assignment
#################

#Use C++11
set(CMAKE_CXX_FLAGS "--std=c++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

#################
#Find CUDA
################
find_package(CUDA REQUIRED)


#################
#Find OpenMPI
################
find_package(MPI)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS})
include_directories(${MPI_INCLUDE_PATH})


############
#Find HDF5 Package for file IO
############

find_package(HDF5 COMPONENTS HL)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_LIBRARIES})

############
#Only if CUDA is found configuration works!
###########

if(${CUDA_FOUND})
  message("-- Using CUDA")
  include_directories("${PROJECT_SOURCE_DIR}/include")
  file(GLOB SRC_FILES "main.cu")
  cuda_add_executable(1dim ${SRC_FILES})
  target_link_libraries(1dim ${LIBS} )

elseif(NOT ${CUDA_FOUND})

  message(FATAL_ERROR "CUDA NOT FOUND BUT IS REQUIRED")

 endif(${CUDA_FOUND})


