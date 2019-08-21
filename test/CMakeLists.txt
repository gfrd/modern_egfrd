# Disable warning 4127 (Conditional Expression is Constant) artifact of TinyTest
if (WIN32)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4127")
endif ()

# Add sub-folders
add_subdirectory(TestGFRD)
add_subdirectory(TestGreensFunctions)
add_subdirectory(TestCCD)


