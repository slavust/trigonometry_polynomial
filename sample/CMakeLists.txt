cmake_minimum_required(VERSION 3.7)

project(sample)
set(CMAKE_CXX_STANDARD 20)
set(SOURCES sample.cpp)
include_directories(../include)
add_executable(sample ${SOURCES})

