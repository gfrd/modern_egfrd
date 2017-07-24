file(GLOB headers *.hpp)
file(GLOB sources *.cpp)

add_executable(TestGFRD ${sources} ${headers})

target_include_directories (TestGFRD PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(TestGFRD PUBLIC eGFRD GreensFunctions Logger ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
