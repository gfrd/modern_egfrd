set (SRC RunGfrd.cpp)
set (HDR Custom.hpp Equilibrium.hpp getoptions.hpp mapK.hpp PowerSpectrum.hpp Resume.hpp Reversible.hpp RunGfrd.hpp Simulation.hpp Simulation3P.hpp ParseIni/CopyNumbersSection.hpp ParseIni/ParserExceptions.hpp ParseIni/ParticlesSection.hpp ParseIni/ReactionRuleSection.hpp ParseIni/SectionBase.hpp ParseIni/SimulatorSection.hpp ParseIni/SimulatorSettings.hpp ParseIni/SpeciesTypeSection.hpp ParseIni/WorldSection.hpp )

add_executable(RunGfrd ${SRC} ${HDR})

if (WIN32)
target_compile_definitions(RunGfrd PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

target_link_libraries(RunGfrd PRIVATE eGFRD)

target_compile_features (RunGfrd PRIVATE cxx_noexcept cxx_deleted_functions cxx_nullptr cxx_auto_type cxx_range_for cxx_override cxx_alias_templates)
