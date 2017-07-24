file(GLOB_RECURSE headers *.hpp)
file(GLOB sources *.cpp)

if (WIN32)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

add_executable(RunGfrd ${sources} ${headers})

target_include_directories (RunGfrd PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(RunGfrd PUBLIC eGFRD GreensFunctions Logger ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})
