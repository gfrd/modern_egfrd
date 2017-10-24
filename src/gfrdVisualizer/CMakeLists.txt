find_package(OpenGL)
find_package(GLUT)

if (OPENGL_FOUND AND GLUT_FOUND)

set (SRC ExternalSim.cpp gfrdRenderer.cpp gfrdVisualizer.cpp ScreenShot.cpp)
set (HDR CameraController.hpp ExternalSim.hpp gfrdRenderer.hpp ScreenShot.hpp)

add_executable(gfrdVisualizer ${SRC} ${HDR})

target_link_libraries(gfrdVisualizer PRIVATE eGFRD ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} GLUT::GLUT)

target_compile_features (gfrdVisualizer PRIVATE cxx_noexcept cxx_deleted_functions cxx_nullptr cxx_auto_type cxx_range_for cxx_override)

if (WIN32)
    target_link_libraries(gfrdVisualizer PUBLIC gdiplus)
	target_compile_definitions(gfrdVisualizer PRIVATE _CRT_SECURE_NO_WARNINGS)
	target_sources(gfrdVisualizer PRIVATE resource.rc)
endif()

else()

message (WARNING, "target gfrdVisualizer skipped")

endif()
