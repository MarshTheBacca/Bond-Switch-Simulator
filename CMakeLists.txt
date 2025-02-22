# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.20)

# Set the project name and the language
project(Bond-Switch-Simulator LANGUAGES CXX)

# Set the C++ standard to C++23
set(CMAKE_CXX_STANDARD 23)

# Require the C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make directory for executable to be saved
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Add include directories
add_subdirectory(src)