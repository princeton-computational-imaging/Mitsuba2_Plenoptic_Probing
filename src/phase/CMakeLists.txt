set(MTS_PLUGIN_PREFIX "phasefunctions")

add_plugin(hg hg.cpp)
add_plugin(isotropic isotropic.cpp)
add_plugin(small_mie small_mie.cpp)

# Register the test directory
add_tests(${CMAKE_CURRENT_SOURCE_DIR}/tests)
