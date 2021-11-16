cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(chain)

# Build target module
add_executable(chain main.cpp flow.h)
# Link Libraries
target_link_libraries(chain tbb)
# Inlcude
target_include_directories(chain PRIVATE .)
# Link s0s
target_include_directories(chain PRIVATE ../../../include)
target_include_directories(chain PRIVATE ../../../../s0s/include)
target_include_directories(chain PRIVATE ../../../../m0sh/include)
target_include_directories(chain PRIVATE ../../../../p0l/include)
target_include_directories(chain PRIVATE ../../../../../thirdparty/eigen)
