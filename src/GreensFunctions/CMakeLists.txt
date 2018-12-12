set (SRC CylindricalBesselGenerator.cpp freeFunctions.cpp GreensFunction1DAbsAbs.cpp GreensFunction1DAbsSinkAbs.cpp GreensFunction1DRadAbs.cpp GreensFunction2DAbsSym.cpp GreensFunction2DRadAbs.cpp GreensFunction3D.cpp GreensFunction3DAbs.cpp GreensFunction3DAbsSym.cpp GreensFunction3DRadAbs.cpp GreensFunction3DRadInf.cpp GreensFunction3DSym.cpp helperFunctionsGf.cpp SphericalBesselGenerator.cpp )
set (HDR CylindricalBesselGenerator.hpp freeFunctions.hpp GreensFunction1DAbsAbs.hpp GreensFunction1DAbsSinkAbs.hpp GreensFunction1DRadAbs.hpp GreensFunction2DAbsSym.hpp GreensFunction2DRadAbs.hpp GreensFunction3D.hpp GreensFunction3DAbs.hpp GreensFunction3DAbsSym.hpp GreensFunction3DRadAbs.hpp GreensFunction3DRadInf.hpp GreensFunction3DSym.hpp helperFunctionsGf.hpp SphericalBesselGenerator.hpp DefsGf.hpp GreensFunction.hpp PairGreensFunction.hpp )

add_library(GreensFunctions SHARED ${SRC} ${HDR})

target_include_directories (GreensFunctions INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} PRIVATE ${GFRD_COMPAT_DIR})

if (WIN32)
target_compile_options(GreensFunctions PRIVATE /wd4251)
endif()

target_link_libraries (GreensFunctions PUBLIC Logger GSL::gsl)

target_compile_features (GreensFunctions PUBLIC cxx_noexcept cxx_deleted_functions cxx_override cxx_alias_templates PRIVATE cxx_nullptr cxx_auto_type cxx_lambdas)
