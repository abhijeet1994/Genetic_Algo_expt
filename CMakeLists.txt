cmake_minimum_required(VERSION 3.14)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(EFL)

# Find packages
find_package(Pagmo REQUIRED)

#Bring the headers into the project
include_directories(include ${PAGMO_INCLUDE_DIRS})

#Generate Shared library from sources
set (FoodData src/getFoodData.cpp)
add_library(getFoodData SHARED ${FoodData})

set (GAlgo src/genetic_algo.cpp)
add_library(geneticAlgo SHARED ${GAlgo})

#Main executable
add_executable(EFLtest src/main.cpp)
target_link_libraries(EFLtest PUBLIC getFoodData geneticAlgo )

#Pagmo test executable
add_executable(PAGMOtest src/tutorial1.cpp)
target_link_libraries(PAGMOtest getFoodData geneticAlgo Pagmo::pagmo)