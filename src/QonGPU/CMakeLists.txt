######
#Right CMake Version
######

cmake_minimum_required(VERSION 3.0)



######
set(CMAKE_COLOR_MAKEFILE ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpedantic")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnarrowing")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Waddress")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Warray-bounds")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wchar-subscripts")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcomment")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wformat")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmain")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmissing-braces")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wparentheses")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsequence-point")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsign-compare")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-aliasing")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wstrict-overflow=1")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wswitch")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wtrigraphs")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wuninitialized")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-function")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-label")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-value")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wvolatile-register-var")

#####
#Set testing options
####
option(TOTEST "Enable test build!" OFF)

######
#Specify Project Name
######

project(qsolve)

######
#Set Compiler Flags
#####

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


######
#Find Packages
#####

#Find Cuda
find_package(CUDA REQUIRED)
#set(LIBS ${LIBS} ${CUDA_LIBRARIES})
#set(LIBS ${LIBS} ${CUDA_CUBLAS_LIBRARIES})
set(LIBS ${LIBS} ${CUDA_cusparse_LIBRARY})
set(CUDA_CXX_FLAGS "${CUDA_CXX_FLAGS} -arch=sm_35")
message("${CUDA_cusparse_LIBRARY}")
#message("${CUDA_LIBRARIES}")
#Find MPI
#find_package(MPI REQUIRED)
#set(LIBS "${LIBS} ${MPI_LIBRARIES}")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS})

#Find HDF5 Package for file IO
FIND_PACKAGE(HDF5 COMPONENTS HL )
include_directories(${HDF5_INCLUDE_DIRS})

set(LIBS ${LIBS} ${HDF5_LIBRARIES})

#find Package Boost
find_package(Boost COMPONENTS unit_test_framework )
include_directories(${BOOST_INCLUDE_DIR})
set(LIBS ${LIBS} ${BOOST_LIBRARIES})

# find openMP
#find_package(OpenMP REQUIRED)
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

ADD_DEFINITIONS(-D__STRICT_ANSI__)
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

#################
#  Enale Testing
#################
if(TOTEST)
enable_testing()


file(GLOB TEST_SRCS_DYN ${CMAKE_CURRENT_SOURCE_DIR}/include/dynamic/tests/*.cu)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_MWAITXINTRIN_H_INCLUDED")

# for each test file set an Test
foreach(testSrc ${TEST_SRCS_DYN})
    get_filename_component(testName ${testSrc} NAME_WE)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    cuda_add_executable( ${testName} ${testSrc})
    target_link_libraries(${testName} ${LIBS})

    set_target_properties(${testName} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/testBin)

    add_test(NAME ${testName}
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/testBin
                 COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testBin/${testName} )

endforeach(testSrc)
endif()
