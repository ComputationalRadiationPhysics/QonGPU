cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "./cmake")

find_package(Splash)

add_subdirectory("./src" "./build/Solver")
