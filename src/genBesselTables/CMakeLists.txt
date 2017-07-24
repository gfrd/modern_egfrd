file(GLOB headers *.hpp)
file(GLOB sources *.cpp)

add_executable(genBesselTables ${sources} ${headers})

target_include_directories (genBesselTables PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries (genBesselTables PRIVATE ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
