cmake_minimum_required(VERSION 3.2.0)

project(modern_egfrd CXX)

# GSL is required for almost all targets
find_package(GSL REQUIRED)

# Specify ISO C++14 version 
if (WIN32)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /std:c++latest")
else (UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall -std=c++14")
endif ()

# Specify output folder as 'bin', multi-configuration tools append the configuration name
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

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

# Specify target dependancies
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
