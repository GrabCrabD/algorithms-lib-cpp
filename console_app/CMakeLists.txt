project(s21_console_app)

add_library(${PROJECT_NAME} STATIC s21_console_app.cpp)

target_link_libraries(${PROJECT_NAME} s21_graph_algorithms readline)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(main main.cpp)

target_link_libraries(main ${PROJECT_NAME})

add_custom_target(console_app COMMAND main)

add_dependencies(console_app s21_graph_algorithms ${PROJECT_NAME})
