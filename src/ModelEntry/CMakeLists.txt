set (SRC SectionBase.cpp VariablesSection.cpp SimulatorSettings.cpp)
set (HDR CopyNumbersSection.hpp ParserExceptions.hpp ParticlePositionsSection.hpp ParticlesSection.hpp ProgressSection.hpp ReactionRecordSection.hpp ReactionRuleSection.hpp SectionBase.hpp SimulatorSection.hpp SimulatorSettings.hpp SpeciesTypeSection.hpp StructureTypeSection.hpp VariablesSection.hpp WorldSection.hpp )

add_library(ModelEntry SHARED ${SRC} ${HDR})

if (WIN32)
target_compile_options(ModelEntry PRIVATE /bigobj /wd4251)
endif()

target_include_directories (ModelEntry INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} PRIVATE ${GFRD_COMPAT_DIR})

target_link_libraries (ModelEntry PUBLIC Logger eGFRD)

target_compile_features (ModelEntry PUBLIC cxx_noexcept cxx_deleted_functions cxx_override cxx_alias_templates PRIVATE cxx_nullptr cxx_auto_type cxx_lambdas)
