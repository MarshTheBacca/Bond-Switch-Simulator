# Download and configure Google Test using FetchContent
include(FetchContent)

# Set FetchContent to use the build directory for its operations
set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/_deps)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
)

# Prevent Google Test from overriding compiler/linker options
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Fetch and make Google Test available
FetchContent_MakeAvailable(googletest)

# Add the test executable
add_executable(run_tests
    test_vector_tools.cpp # Add your test files here
)

# Add the include directory to the include paths
target_include_directories(run_tests PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link Google Test and your project to the test executable
target_link_libraries(run_tests
    gtest
    gtest_main
)

# Add tests to CTest
include(GoogleTest)

gtest_discover_tests(run_tests)