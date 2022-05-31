set(MTS_PLUGIN_PREFIX "sensors")

add_plugin(perspective       perspective.cpp)
add_plugin(radiancemeter     radiancemeter.cpp)
add_plugin(thinlens          thinlens.cpp)
add_plugin(irradiancemeter   irradiancemeter.cpp)
#add_plugin(perspective_rdist perspective_rdist.cpp)
#add_plugin(coded_perspective coded_perspective.cpp)
#add_plugin(coded_orthographic coded_orthographic.cpp)

# Register the test directory
add_tests(${CMAKE_CURRENT_SOURCE_DIR}/tests)
