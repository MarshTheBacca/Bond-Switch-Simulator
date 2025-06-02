# Create a new executable
add_executable(bond_switch_simulator)

# Use CMAKE_PREFIX_PATH to search for lammps
find_library(LAMMPS_LIBRARY lammps)

# Get LAMMPS include directory
find_path(LAMMPS_INCLUDE_DIR lammps.h PATHS $ENV{CMAKE_PREFIX_PATH}/include/lammps)

# Check if LAMMPS include directory was found
if(NOT LAMMPS_INCLUDE_DIR)
# Try to set a default value
    message(WARNING "Could not find LAMMPS include directory using CMAKE_PREFIX_PATH. Trying default path.")
    set(LAMMPS_INCLUDE_DIR $ENV{HOME}/.local/include/lammps)
    # Check if the default value exists
    if(NOT EXISTS ${LAMMPS_INCLUDE_DIR})
        message(FATAL_ERROR "LAMMPS include directory not found")
    endif()
endif()


# Check if LAMMPS library was found
if(NOT LAMMPS_LIBRARY)
    # Try to set a default value
    message(WARNING "Could not find LAMMPS library path using CMAKE_PREFIX_PATH. Trying default path.")

    set(LAMMPS_LIBRARY $ENV{HOME}/.local/lib/liblammps.so)
    # Check if the default value exists
    if(NOT EXISTS ${LAMMPS_LIBRARY})
        message(FATAL_ERROR "LAMMPS library not found")
    endif()
endif()

# Print the LAMMPS library path and include directory
message(STATUS "LAMMPS library: ${LAMMPS_LIBRARY}")
message(STATUS "LAMMPS include directory: ${LAMMPS_INCLUDE_DIR}")

message(STATUS "Finding spdlog and MPI...")
find_package(spdlog REQUIRED)
find_package(MPI REQUIRED)

# Add all source files to the executable
file(GLOB SOURCES "*.cpp")
target_sources(bond_switch_simulator PRIVATE ${SOURCES})

# Add include directories
target_include_directories(bond_switch_simulator
    BEFORE PRIVATE ${CMAKE_SOURCE_DIR}/include
    PRIVATE ${LAMMPS_INCLUDE_DIR}
    PRIVATE ${MPI_CXX_INCLUDE_DIRS}
)
# Link the LAMMPS, spdlog, and MPI libraries
target_link_libraries(bond_switch_simulator
    PRIVATE ${LAMMPS_LIBRARY}
    PRIVATE spdlog::spdlog_header_only
    PRIVATE MPI::MPI_CXX
)