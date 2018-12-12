set (SRC RunGfrd.cpp)
set (HDR ../Common/getoptions.hpp SimModel.hpp SimResume.hpp SimCustom.hpp Simulation.hpp  )

add_executable(RunGfrd ${SRC} ${HDR})

target_link_libraries(RunGfrd PRIVATE eGFRD ModelEntry)

target_compile_features (RunGfrd PRIVATE cxx_noexcept cxx_deleted_functions cxx_nullptr cxx_auto_type cxx_range_for cxx_override cxx_alias_templates)

if (WIN32)
	target_compile_definitions(RunGfrd PRIVATE _CRT_SECURE_NO_WARNINGS)
	target_compile_options(RunGfrd PRIVATE /wd4251)

	configure_file(Resources/resource.rc.in resource.rc @ONLY)
	target_include_directories(RunGfrd PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Resources )
	target_sources(RunGfrd PRIVATE Resources/resource.h ${CMAKE_CURRENT_BINARY_DIR}/resource.rc)
endif()
