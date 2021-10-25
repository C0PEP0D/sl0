cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(shear)

# Used for Vim
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )
IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
        )
ENDIF()

# Build target module
add_executable(shear main.cpp)
target_include_directories(shear PRIVATE .)

target_include_directories(shear PRIVATE ../../../sa0/include)
target_include_directories(shear PRIVATE ../../../s0s/include)
target_include_directories(shear PRIVATE ../../../sl0/include)
target_include_directories(shear PRIVATE ../../../fl0w/include)
target_include_directories(shear PRIVATE ../../../../thirdparty/eigen)
