file(GLOB headers *.hpp)
file(GLOB sources *.cpp)

if (WIN32)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()
add_definitions(-DLOGGER_EXPORTS)

add_library(Logger SHARED ${sources} ${headers})

target_include_directories (Logger INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})