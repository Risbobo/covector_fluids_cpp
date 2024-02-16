cmake_minimum_required(VERSION 3.25)
project(covector_fluids_cpp)

set(CMAKE_CXX_STANDARD 20)

add_executable(covector_fluids_cpp main.cpp StGrid.cpp StGrid.h FlowCpp.cpp FlowCpp.h
        StGridv2.cpp
        StGridv2.h
        FlowCppv2.cpp
        FlowCppv2.h)
