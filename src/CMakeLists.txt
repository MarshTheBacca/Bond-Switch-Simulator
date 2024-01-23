set(SOURCE_FILES
    lammps_c_interface.cpp
    linked_network.cpp
    main.cpp
    monte_carlo.cpp
    network.cpp
    node.cpp
    opt.tpp
    pot.cpp
    pot1d.cpp
    pot2d.cpp
    pot3d.cpp
    vecf.tpp
    vecr.tpp
)

link_directories(${LAMMPS_LIBRARY_DIRS})
add_executable(netmc.x ${SOURCE_FILES})
target_include_directories(netmc.x BEFORE PRIVATE ${LAMMPS_INCLUDE_DIRS})
target_link_libraries(netmc.x PRIVATE ${LAMMPS_LIBRARIES} ${OpenMP_C_FLAGS} ${OpenMP_CXX_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})

message(STATUS "LAMMPS include directories: ${LAMMPS_INCLUDE_DIRS}")