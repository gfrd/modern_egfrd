if (APPLE)
	# On Mac OS-X we need FreeGLUT and not the default system GLUT which is found by find_package(GLUT).
	# furthermore FreeGLUT (from homebrew) links to the XQuartz version of OpenGL, so we also need to check that.

	# specify your Homebrew prefix (brew --prefix)
	set (_GLUT_ROOT_PATH    "/usr/local/"  )

	# find FreeGLUT
	find_path (_GLUT_INCLUDE_DIR GL/freeglut.h PATHS ${_GLUT_ROOT_PATH}/include)
	find_library (_GLUT_LIBRARY libglut.dylib PATHS ${_GLUT_ROOT_PATH}/lib NO_DEFAULT_PATH )
	
	if (NOT _GLUT_INCLUDE_DIR OR NOT _GLUT_LIBRARY)
		message (WARNING, "FreeGLUT was not found!")
	else()
		message("Found FreeGLUT (include: ${_GLUT_INCLUDE_DIR} and lib: ${_GLUT_LIBRARY} )")

		set (GLUT_FOUND TRUE)
		add_library(GLUT::GLUT UNKNOWN IMPORTED)
		set_target_properties(GLUT::GLUT PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_GLUT_INCLUDE_DIR}")
		set_target_properties(GLUT::GLUT PROPERTIES IMPORTED_LOCATION "${_GLUT_LIBRARY}")

		# If FreeGLUT was installed with brew, it uses the XQuartz window manager
		# otherwise it most likely links to the native window manager.
		# Here you have to choose! No automated test available (maybe check pkg-config ??).
		if (NOT TRUE)
			# find native version of OpenGL
			find_package(OpenGL)
		else()
			# find XQuartz/OpenGL (lines copied from freeglut.rb)
			find_path(_GLX_INCLUDE_DIR GL/glx.h PATHS /opt/X11/include /usr/X11/include /usr/X11R6/include)
			find_library(_OPENGL_gl_LIBRARY GL PATHS /opt/X11/lib /usr/X11/lib /usr/X11R6/lib)
			find_library(_OPENGL_glu_LIBRARY GLU PATHS /opt/X11/lib /usr/X11/lib /usr/X11R6/lib)
			
			if (NOT _GLX_INCLUDE_DIR OR NOT _OPENGL_gl_LIBRARY)
				message (WARNING, "OpenGL(XQuartz) was not found!")
			else()
				message("Found OpenGL (include: ${_GLX_INCLUDE_DIR} and libs: ${_OPENGL_gl_LIBRARY} and ${_OPENGL_glu_LIBRARY} )")

				set (OPENGL_FOUND TRUE)
				list(APPEND OPENGL_LIBRARIES ${_OPENGL_gl_LIBRARY} ${_OPENGL_glu_LIBRARY})
				set_target_properties(GLUT::GLUT PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_GLX_INCLUDE_DIR}")
			endif()
		endif()
	endif()
	
else()
	find_package(OpenGL)
	find_package(GLUT)
endif()


if (OPENGL_FOUND AND GLUT_FOUND)

set (SRC ExternalSim.cpp gfrdRenderer.cpp gfrdVisualizer.cpp ScreenShot.cpp)
set (HDR ../Common/getoptions.hpp CameraController.hpp ExternalSim.hpp gfrdRenderer.hpp ScreenShot.hpp)

add_executable(gfrdVisualizer ${SRC} ${HDR})

target_link_libraries(gfrdVisualizer PRIVATE eGFRD ModelEntry ${OPENGL_LIBRARIES} GLUT::GLUT)

target_compile_features (gfrdVisualizer PRIVATE cxx_noexcept cxx_deleted_functions cxx_nullptr cxx_auto_type cxx_range_for cxx_override)

if (WIN32)
	target_link_libraries(gfrdVisualizer PUBLIC gdiplus)
	target_compile_definitions(gfrdVisualizer PRIVATE _CRT_SECURE_NO_WARNINGS)
	target_compile_options(gfrdVisualizer PRIVATE /wd4251)
	configure_file(Resources/resource.rc.in resource.rc @ONLY)
	target_include_directories(gfrdVisualizer PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Resources )
	target_sources(gfrdVisualizer PRIVATE Resources/resource.h ${CMAKE_CURRENT_BINARY_DIR}/resource.rc)
elseif (APPLE)
#	target_compile_options(gfrdVisualizer -Wno-deprecated-declarations)
endif()


else()

message (WARNING, "target gfrdVisualizer skipped")

endif()
