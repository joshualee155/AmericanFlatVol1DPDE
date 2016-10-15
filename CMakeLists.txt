cmake_minimum_required(VERSION 3.6)
project(AmericanFlatVol1DPDE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp clsAmericanFlatVol1DPDE.cpp)
add_executable(AmericanFlatVol1DPDE ${SOURCE_FILES})