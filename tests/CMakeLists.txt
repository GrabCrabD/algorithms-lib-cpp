add_subdirectory(graph_algorithms_tests)
add_subdirectory(graph_tests)

## Находим gcovr
find_program(GCOVR_EXECUTABLE gcovr)

if(COVERAGE AND GCOVR_EXECUTABLE)
    message(STATUS "Измерение покрытия кода тестами включено")

    ## Создаем команду для запуска анализа покрытия тестами
  add_custom_target(coverage
        COMMAND 
        ${GCOVR_EXECUTABLE}
        --root=${CMAKE_SOURCE_DIR}/source # указываем папки с анализируемыми исходниками
        --object-directory=${CMAKE_BINARY_DIR} # указываем папку с бинарниками
        --html --html-details -o ${CMAKE_SOURCE_DIR}/report/report.html # просим создать отчет в HTML
        # говорим, что данная цель зависит от команды tests
        DEPENDS graph_tests graph_algorithms_tests
   )
elseif(MYLIB_COVERAGE AND NOT GCOVR_EXECUTABLE)
    set(MYLIB_COVERAGE OFF)
    message(WARNING "Для замеров покрытия кода тестами требуется программа gcovr")
endif()
