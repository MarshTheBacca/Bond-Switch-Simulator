# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.20)

# Set the project name and the language
project(Bond-Switch-Simulator LANGUAGES CXX)

# Set the C++ standard to C++23
set(CMAKE_CXX_STANDARD 23)

# Require the C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Enable color diagnostics for Clang
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-fcolor-diagnostics)
endif()

# Enable color output for Ninja
if (CMAKE_GENERATOR STREQUAL "Ninja")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
endif()

# Make directory for executable to be saved
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Add include directories
add_subdirectory(src)

# Add the tests directory
enable_testing()
add_subdirectory(tests)