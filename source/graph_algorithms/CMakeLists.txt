project(s21_graph_algorithms)

add_library(${PROJECT_NAME} STATIC s21_graph_algorithms.cpp)

target_link_libraries(${PROJECT_NAME} s21_graph s21_queue s21_stack)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
