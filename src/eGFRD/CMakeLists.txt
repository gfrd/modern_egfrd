find_package(ccd CONFIG REQUIRED)

file(GLOB SRC EGFRDSimulator.cpp Matrix4.cpp Multi.cpp Pair.cpp Persistence.cpp PlanarSurface.cpp scaling.cpp Single.cpp Vector2.cpp Vector3.cpp)
file(GLOB HDR BDPropagator.hpp BDSimulator.hpp Box.hpp ConnectivityContainer.hpp convert.hpp CopyNumbers.hpp CuboidalRegion.hpp Cylinder.hpp CylindricalSurface.hpp DefsEgfrd.hpp Disk.hpp DiskSurface.hpp Domain.hpp DomainID.hpp DynamicPriorityQueue.hpp EGFRDSimulator.hpp Event.hpp EventID.hpp EventScheduler.hpp exceptions.hpp genericIterator.hpp GreenFunctionHelpers.hpp helperFunctions.hpp Identifier.hpp InteractionRule.hpp makeString.hpp Matrix4.hpp MatrixSpace.hpp Model.hpp Multi.hpp Pair.hpp Particle.hpp ParticleContainer.hpp ParticleContainerImpl.hpp ParticleID.hpp ParticleSimulator.hpp Persistence.hpp PlanarSurface.hpp Plane.hpp Progress.hpp randomNumberGenerator.hpp Reactants.hpp ReactionRecorder.hpp ReactionRule.hpp ReactionRuleCollection.hpp ReactionRuleID.hpp Region.hpp scaling.hpp SerialIDGenerator.hpp Shell.hpp ShellCreateUtils.hpp ShellID.hpp Single.hpp SpeciesType.hpp SpeciesTypeID.hpp Sphere.hpp SphericalSurface.hpp Structure.hpp StructureContainer.hpp StructureID.hpp StructureType.hpp StructureTypeID.hpp Surface.hpp Transaction.hpp Vector2.hpp Vector3.hpp World.hpp)

add_library(eGFRD SHARED ${SRC} ${HDR})

target_include_directories (eGFRD INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries (eGFRD PUBLIC GreensFunctions ccd)

target_compile_features (eGFRD PUBLIC cxx_strong_enums cxx_noexcept cxx_deleted_functions cxx_nullptr cxx_auto_type cxx_range_for cxx_rvalue_references cxx_override cxx_lambdas cxx_alias_templates )

if(!APPLE)
    set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined" )
endif()