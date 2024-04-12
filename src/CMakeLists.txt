find_package(PkgConfig REQUIRED)
pkg_check_modules(LAMMPS REQUIRED IMPORTED_TARGET liblammps)

# Find the OpenMP package
find_package(OpenMP REQUIRED)

find_package(spdlog REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(bond_switch_simulator.exe)
target_sources(bond_switch_simulator.exe PRIVATE
    lammps_object.cpp
    linked_network.cpp
    main.cpp
    metropolis.cpp
    network.cpp
    node.cpp
    input_data.cpp
    output_file.cpp
    vector_tools.cpp
)
target_include_directories(bond_switch_simulator.exe PUBLIC ${LAMMPS_INCLUDE_DIRS}/lammps)
target_include_directories(bond_switch_simulator.exe BEFORE PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link against the LAMMPS and OpenMP libraries
target_link_libraries(bond_switch_simulator.exe PRIVATE PkgConfig::LAMMPS OpenMP::OpenMP_CXX spdlog::spdlog_header_only)

message(STATUS "LAMMPS include directories: ${LAMMPS_INCLUDE_DIRS}")

