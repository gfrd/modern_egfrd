set (SRC genBesselTables.cpp)
set (HDR "")

add_executable(genBesselTables ${SRC} ${HDR})

target_link_libraries (genBesselTables PRIVATE GSL::gsl)

target_compile_features (genBesselTables PRIVATE cxx_alias_templates )

