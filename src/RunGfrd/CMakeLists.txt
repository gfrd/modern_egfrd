set (SRC RunGfrd.cpp ParseIni/SectionBase.cpp ParseIni/VariablesSection.cpp)
set (HDR getoptions.hpp SimModel.hpp SimResume.hpp Simulation.hpp ParseIni/CopyNumbersSection.hpp ParseIni/ParserExceptions.hpp ParseIni/ParticlePositionsSection.hpp ParseIni/ParticlesSection.hpp ParseIni/ReactionRecordSection.hpp ParseIni/ReactionRuleSection.hpp ParseIni/SectionBase.hpp ParseIni/SimulatorSection.hpp ParseIni/SimulatorSettings.hpp ParseIni/SpeciesTypeSection.hpp ParseIni/WorldSection.hpp ParseIni/ProgressSection.hpp )

add_executable(RunGfrd ${SRC} ${HDR})

if (WIN32)
target_compile_definitions(RunGfrd PRIVATE _CRT_SECURE_NO_WARNINGS _SCL_SECURE_NO_WARNINGS)
target_compile_options(RunGfrd PRIVATE /bigobj)
endif()

target_link_libraries(RunGfrd PRIVATE eGFRD)

target_compile_features (RunGfrd PRIVATE cxx_noexcept cxx_deleted_functions cxx_nullptr cxx_auto_type cxx_range_for cxx_override cxx_alias_templates)
