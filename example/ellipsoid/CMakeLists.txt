cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(ellipsoid)

# Build target module
add_executable(spheroid main.cpp flow.h)
target_include_directories(spheroid PRIVATE .)

# Link s0s
target_include_directories(spheroid PRIVATE ../../include)
target_include_directories(spheroid PRIVATE ../../../s0s/include)
target_include_directories(spheroid PRIVATE ../../../../thirdparty/eigen)
