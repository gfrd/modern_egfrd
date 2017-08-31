set (SRC CylindricalBesselGenerator_test.cpp Factorial_test.cpp GreensFunction3D_test.cpp GreensFunction3DAbs_test.cpp GreensFunction3DAbsSym_test.cpp GreensFunction3DRadAbs_test.cpp GreensFunction3DRadInf_test.cpp GslSumLevinU_test.cpp SphericalBesselGenerator_test.cpp tinytest.cpp ../common/tinytest_specific.cpp)
set (HDR CylindricalBesselGenerator_test.hpp Factorial_test.hpp GreensFunction3D_test.hpp GreensFunction3DAbs_test.hpp GreensFunction3DAbsSym_test.hpp GreensFunction3DRadAbs_test.hpp GreensFunction3DRadInf_test.hpp GslSumLevinU_test.hpp SphericalBesselGenerator_test.hpp ../common/tinytest.h)

add_executable(TestGreensFunctions ${SRC} ${HDR})

target_link_libraries(TestGreensFunctions PRIVATE GreensFunctions)

target_compile_features (TestGreensFunctions PRIVATE cxx_auto_type)
