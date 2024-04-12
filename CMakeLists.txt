cmake_minimum_required(VERSION 3.10)

project(NetMCMarsh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force C++17 or better")
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add an option for the executable export path
set(EXEC_OUTPUT_PATH ${CMAKE_BINARY_DIR} CACHE PATH "Path to put the executable")

# Make directory for executable to be saved
file(MAKE_DIRECTORY ${EXEC_OUTPUT_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXEC_OUTPUT_PATH})

# Add include directories
add_subdirectory(src)