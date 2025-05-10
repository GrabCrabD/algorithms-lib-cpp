include(GoogleTest)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR})

add_executable(graph_tests graph_tests.cpp)

target_link_libraries(graph_tests s21_graph GTest::gtest_main)
gtest_discover_tests(graph_tests)

add_custom_target(graph_testing COMMAND graph_tests)
add_dependencies(graph_testing s21_graph)
