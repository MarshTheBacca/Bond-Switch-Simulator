find_package(PkgConfig REQUIRED)
pkg_check_modules(LAMMPS REQUIRED IMPORTED_TARGET liblammps)

# Find the OpenMP package
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(netmc.x)
target_sources(netmc.x PRIVATE
    lammps_c_interface.cpp
    linked_network.cpp
    main.cpp
    monte_carlo.cpp
    network.cpp
    node.cpp
    pot.cpp
    pot1d.cpp
    pot2d.cpp
    pot3d.cpp
)

target_include_directories(netmc.x BEFORE PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link against the LAMMPS and OpenMP libraries
target_link_libraries(netmc.x PRIVATE PkgConfig::LAMMPS OpenMP::OpenMP_CXX)

message(STATUS "LAMMPS include directories: ${LAMMPS_INCLUDE_DIRS}")