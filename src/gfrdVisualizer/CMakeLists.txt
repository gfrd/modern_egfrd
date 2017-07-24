find_package(OpenGL)
find_package(GLUT)

if (OPENGL_FOUND AND GLUT_FOUND)

file(GLOB headers *.hpp)
file(GLOB sources *.cpp)

if (WIN32)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    enable_language(RC)
    set(sources Resource.rc ${sources})
endif()

add_executable(gfrdVisualizer ${sources} ${headers})

target_include_directories (gfrdVisualizer PRIVATE ${GSL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIR})
target_link_libraries(gfrdVisualizer PUBLIC eGFRD GreensFunctions Logger ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})

if (WIN32)
    target_link_libraries(gfrdVisualizer PUBLIC gdiplus)
endif()

else()

message (WARNING, "target gfrdVisualizer skipped")

endif()
