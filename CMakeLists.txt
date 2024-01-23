cmake_minimum_required(VERSION 3.10)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

project(NetMCMarsh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force C++11 or better")
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Make directory for executable to be saved, this is where you can put netmc.inpt and such
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)

# Find OpenMP and LAMMPS
# PkgConfig is used to find LAMMPS no matter where it is installed by the user
find_package(PkgConfig REQUIRED)
find_package(OpenMP REQUIRED)
pkg_check_modules(LAMMPS REQUIRED liblammps)

# Add include directories
add_subdirectory(src)

add_custom_target(clean-cmake-files
   COMMAND ${CMAKE_COMMAND} -P clean-all.cmake
)

# Create a clean-all.cmake script
file(WRITE clean-all.cmake "
file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/CMakeCache.txt)
file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/CMakeFiles)
file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/*.cmake)
file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/*.o)
file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/*.ghc)
")

