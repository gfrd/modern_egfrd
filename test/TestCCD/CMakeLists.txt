find_package(ccd CONFIG REQUIRED)

set (SRC LibCCD_test.cpp tinytest.cpp ../common/tinytest_specific.cpp)
set (HDR LibCCD_test.hpp ../common/tinytest.h)

add_executable(TestCCD ${SRC} ${HDR})

target_link_libraries(TestCCD PRIVATE ccd eGFRD)

target_compile_features (TestCCD PRIVATE cxx_auto_type)
