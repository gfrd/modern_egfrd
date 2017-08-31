set (SRC SeedGenerator.cpp)
set (HDR "")

add_executable(SeedGenerator ${SRC} ${HDR})

target_compile_features (SeedGenerator PRIVATE cxx_range_for)
