cmake_minimum_required(VERSION 3.7)

project(benchmark)
set(CMAKE_CXX_STANDARD 20)
set(SOURCES generate_stats.cpp)
include_directories(../include)
add_executable(benchmark ${SOURCES})

