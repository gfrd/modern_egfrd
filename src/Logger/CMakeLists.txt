set (SRC Logger.cpp)
set (HDR Logger.hpp bitmask_operators.hpp)

add_library(Logger SHARED ${SRC} ${HDR})

if (WIN32)
target_compile_definitions(Logger PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

target_include_directories (Logger INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} PUBLIC ${GFRD_COMPAT_DIR})

target_compile_features (Logger PUBLIC cxx_strong_enums cxx_noexcept cxx_deleted_functions cxx_nullptr cxx_rvalue_references PRIVATE cxx_auto_type )
