file(GLOB headers *.hpp)
file(GLOB sources *.cpp)

add_executable(SeedGenerator ${sources} ${headers})
