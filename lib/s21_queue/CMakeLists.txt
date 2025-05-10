project(s21_queue)

add_library(${PROJECT_NAME} INTERFACE)

target_link_libraries(
    ${PROJECT_NAME}
    INTERFACE
        s21_list)

target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
