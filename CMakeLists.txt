cmake_minimum_required(VERSION 3.10)

project(NetMCMarsh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force C++17 or better")
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Make directory for executable to be saved, this is where you can put netmc.inpt and such
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)

# Add include directories
add_subdirectory(src)
