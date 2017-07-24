file(GLOB headers *.hpp)
file(GLOB sources *.cpp)

add_executable(TestGreensFunctions ${sources} ${headers})

target_include_directories (TestGreensFunctions PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(TestGreensFunctions PUBLIC eGFRD GreensFunctions Logger ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
