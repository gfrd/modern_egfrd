set (SRC BDSimulator_test.cpp Box_test.cpp Cylinder_test.cpp Disk_test.cpp EventScheduler_test.cpp GenericIterator_test.cpp Identifier_test.cpp MatrixSpace_test.cpp Model_test.cpp ParticleContainer_test.cpp ParticleSimulator_test.cpp Persistence_test.cpp Plane_test.cpp Single_test.cpp Sphere_test.cpp TestSimple.cpp tinytest.cpp Vector2_test.cpp Vector3_test.cpp World_test.cpp tinytest.cpp ../common/tinytest_specific.cpp)
set (HDR BDSimulator_test.hpp Box_test.hpp Cylinder_test.hpp Disk_test.hpp EventScheduler_test.hpp GenericIterator_test.hpp Identifier_test.hpp MatrixSpace_test.hpp Model_test.hpp ParticleContainer_test.hpp ParticleSimulator_test.hpp Persistence_test.hpp Plane_test.hpp Single_test.hpp Sphere_test.hpp TestSimple.hpp Vector2_test.hpp Vector3_test.hpp World_test.hpp ../common/tinytest.h)

add_executable(TestGFRD ${SRC} ${HDR})

target_link_libraries(TestGFRD PRIVATE eGFRD)

target_compile_features (TestGFRD PRIVATE cxx_nullptr cxx_auto_type cxx_range_for cxx_rvalue_references cxx_alias_templates)
