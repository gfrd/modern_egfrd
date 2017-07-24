file(GLOB headers *.hpp)
file(GLOB sources *.cpp)

add_definitions(-DEGFRD_EXPORTS)

add_library(eGFRD SHARED ${sources} ${headers})

target_include_directories (eGFRD PRIVATE ${GSL_INCLUDE_DIRS})
target_include_directories (eGFRD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (eGFRD PUBLIC GreensFunctions Logger ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
