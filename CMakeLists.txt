cmake_minimum_required(VERSION 3.16)
project(graphs VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

#Задааем переменную для покрытия тестами и выключаем тестирование по дефолту
option(COVERAGE "Покрытие тестами" OFF)

## Проверяем если включено покрытие тестами(gcov), то добавляем флаг компиляции
if(COVERAGE)
add_compile_options("--coverage")
add_link_options("--coverage")
endif()

#Задаем флаги компиляции
if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
add_compile_options(
    -Wall
    -Wextra
    -Wno-deprecated #флаг игнорирующий ошибку устаревших функций
    -Wpedantic
)
endif()

add_subdirectory(lib) #собираем либы стэка и очереди
add_subdirectory(source/graph) #собираем либу парсинга графа из файла и экспорта в файл
add_subdirectory(source/graph_algorithms) #собираем либу инкапсулирующую методы работы с графами
add_subdirectory(console_app) #собираем исполняемый файл консольного приложения для работы с графами

add_subdirectory(tests EXCLUDE_FROM_ALL) #собираем и запускаем тесты
