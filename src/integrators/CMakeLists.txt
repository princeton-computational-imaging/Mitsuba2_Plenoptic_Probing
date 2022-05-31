set(MTS_PLUGIN_PREFIX "integrators")

add_plugin(depth   depth.cpp)
add_plugin(direct  direct.cpp)
add_plugin(path    path.cpp)
add_plugin(spatial_probing  spatial_probing.cpp)
add_plugin(aov     aov.cpp)
add_plugin(stokes  stokes.cpp)
add_plugin(moment  moment.cpp)
add_plugin(volpath  volpath.cpp)
add_plugin(time_gate  time_gate.cpp)

add_plugin(sp_volpath  sp_volpath.cpp)
add_plugin(volpathmis volpathmis.cpp)
add_plugin(volpath_simple volpath_simple.cpp)
add_plugin(null_path null_path.cpp)

# Register the test directory
add_tests(${CMAKE_CURRENT_SOURCE_DIR}/tests)
