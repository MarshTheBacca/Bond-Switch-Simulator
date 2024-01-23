cmake_minimum_required(VERSION 3.10)

project(NetMCMarsh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force C++11 or better")
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Make directory for executable to be saved, this is where you can put netmc.inpt and such
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)

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

