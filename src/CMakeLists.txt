cmake_minimum_required(VERSION 3.16)
project(cfSNV)

if (UNIX)
    set (CMAKE_CXX_FLAGS "-g3 -Wall -O3")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_COMPILER g++)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -c")
set(CMAKE_VERBOSE_MAKEFILE 1)


#include_directories(c)

#add_library(helper STATIC c/_helper.cpp)
#target_include_directories(helper PUBLIC c)

add_executable(cfSNV
        c/_probability.cpp
        c/_helper.h
        c/_helper.cpp c/estimate_TFestimate.cpp c/estimate_TFestimate.h c/_parameter.h c/hotspot.cpp c/hotspot.h c/_probability.h c/_filter.cpp c/_filter.h c/_parameter.cpp)
