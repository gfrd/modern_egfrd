file(GLOB headers *.hpp)
file(GLOB sources *.cpp)

add_definitions(-DGREENSFUNCTIONS_EXPORTS)

add_library(GreensFunctions SHARED ${sources} ${headers})

target_include_directories (GreensFunctions PRIVATE ${GSL_INCLUDE_DIRS})
target_include_directories (GreensFunctions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (GreensFunctions PUBLIC Logger ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
