cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(static)

# Build target module
add_executable(static main.cpp flow.h)
# Link Libraries
target_link_libraries(static tbb)
# Inlcude
target_include_directories(static PRIVATE .)
# Link s0s
target_include_directories(static PRIVATE ../../../include)
target_include_directories(static PRIVATE ../../../../s0s/include)
target_include_directories(static PRIVATE ../../../../../thirdparty/eigen)
