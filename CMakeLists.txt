
cmake_minimum_required(VERSION 3.0.0)

add_library(libsplash "./libs/libsplash/src/include/splash/splash_serial.h")
set_target_properties(libsplash PROPERTIES LINKER_LANGUAGE CXX)
add_subdirectory("./src" "./build")