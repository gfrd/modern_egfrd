cmake_minimum_required(VERSION 3.2.0)

project(modern_egfrd CXX)

# Get the build day (days since 1st Jan 2000)
string(TIMESTAMP _GFRD_YEAR %Y)
string(TIMESTAMP _GFRD_DOY %j)
math (EXPR BUILDDAYNUMBER "365*(${_GFRD_YEAR}-2000) + (${_GFRD_YEAR}-2000+3)/4 + ${_GFRD_DOY}-1" )

# The version number
set (GFRD_VERSION_MAJOR 2)
set (GFRD_VERSION_MINOR 1)
set (GFRD_VERSION_PATCH ${BUILDDAYNUMBER})
string(TIMESTAMP GFRD_VERSION_BUILD)

# GSL is required for almost all targets
find_package(GSL REQUIRED)

# Check C++ stdlib features
include(cmake/CheckFeature.cmake)
check_make_unique()
configure_file (cmake/gfrd_compat.hpp.in gfrd_compat.hpp)
set (GFRD_COMPAT_DIR ${CMAKE_BINARY_DIR})

# Specify output folder as 'bin', multi-configuration tools append the configuration name
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Add sub-folders for src and test
add_subdirectory(src)
add_subdirectory(test)

# Custom command for post-build generation of bessel-lookup-tables
add_custom_command(
   OUTPUT ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/SphericalBesselTable
   COMMAND genBesselTables ARGS -s ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/SphericalBesselTable
   DEPENDS genBesselTables
)
add_custom_command(
   OUTPUT ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CylindricalBesselTable
   COMMAND genBesselTables ARGS -c ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CylindricalBesselTable
   DEPENDS genBesselTables
)
add_custom_target(tables ALL DEPENDS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/SphericalBesselTable
                                     ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CylindricalBesselTable)

# Specify target dependencies
add_dependencies(GreensFunctions Logger tables)
add_dependencies(eGFRD GreensFunctions Logger)
add_dependencies(RunGfrd eGFRD)
if (TARGET gfrdVisualizer)
add_dependencies(gfrdVisualizer eGFRD)
endif()

# Build with different Matrix Space Size?
if (DEFINED MATRIXSPACE)
target_compile_definitions(eGFRD PRIVATE MATRIXSIZE=${MATRIXSPACE})
target_compile_definitions(RunGfrd PRIVATE MATRIXSIZE=${MATRIXSPACE})
target_compile_definitions(TestGFRD PRIVATE MATRIXSIZE=${MATRIXSPACE})
if (TARGET gfrdVisualizer)
target_compile_definitions(gfrdVisualizer PRIVATE MATRIXSIZE=${MATRIXSPACE})
endif()
endif()


# Setup testing
enable_testing()
add_test(NAME test1 COMMAND TestGreensFunctions)
add_test(NAME test2 COMMAND TestGFRD)
add_test(NAME test3 COMMAND TestCCD)
