cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(dynamic)

# Build target module
add_executable(dynamic main.cpp flow.h)
# Link Libraries
target_link_libraries(dynamic tbb)
# Inlcude
target_include_directories(dynamic PRIVATE .)
# Link s0s
target_include_directories(dynamic PRIVATE ../../../include)
target_include_directories(dynamic PRIVATE ../../../../s0s/include)
target_include_directories(dynamic PRIVATE ../../../../../thirdparty/eigen)
