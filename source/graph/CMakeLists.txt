project(s21_graph)

add_library(${PROJECT_NAME} STATIC s21_graph.cpp)

target_link_libraries(${PROJECT_NAME} s21_stack s21_queue)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
