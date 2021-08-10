cmake_minimum_required(VERSION 3.16)
project(cfSNV)

if (UNIX)
    set (CMAKE_CXX_FLAGS "-g3 -Wall")
endif ()

#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/local/opt/boost")
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/local/opt/boost/lib")

set(Boost_DEBUG 1)

find_package( Boost REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_COMPILER g++)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -c")
set(CMAKE_VERBOSE_MAKEFILE 1)
set(D_GLIBCXX_USE_CXX11_ABI,"0")


#include_directories(c)

#add_library(helper STATIC c/_helper.cpp)
#target_include_directories(helper PUBLIC c)

add_executable(tfEstimate
        c/_probability.h c/_probability.cpp
        c/_helper.h c/_helper.cpp
        c/_parameter.h c/_parameter.cpp
        c/hotspot.h c/hotspot.cpp
        c/_filter.h c/_filter.cpp
        c/_jenks.h c/_jenks.cpp
        c/estimate_TFestimate.h c/estimate_TFestimate.cpp)

add_executable(genotype
        c/_probability c/_probability.cpp
        c/_helper.h c/_helper.cpp
        c/_parameter.h c/_parameter.cpp
        c/hotspot.h c/hotspot.cpp
        c/_filter.h c/_filter.cpp
        c/_jenks.h c/_jenks.cpp
        c/genotype.genotype.h c/genotype.genotype.cpp)

add_executable(filterWithPileup
        c/_probability c/_probability.cpp
        c/_helper.h c/_helper.cpp
        c/_parameter.h c/_parameter.cpp
        c/hotspot.h c/hotspot.cpp
        c/_filter.h c/_filter.cpp
        c/_jenks.h c/_jenks.cpp
        c/filter_with_pileup.h c/filter_with_pileup.cpp)

TARGET_LINK_LIBRARIES( tfEstimate LINK_PUBLIC ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES( genotype LINK_PUBLIC ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES( filterWithPileup LINK_PUBLIC ${Boost_LIBRARIES} )