cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(point)

# Build target module
add_executable(point main.cpp flow.h)
target_include_directories(point PRIVATE .)

# Link s0s
target_include_directories(point PRIVATE ../../include)
target_include_directories(point PRIVATE ../../../s0s/include)
target_include_directories(point PRIVATE ../../../../thirdparty/eigen)
