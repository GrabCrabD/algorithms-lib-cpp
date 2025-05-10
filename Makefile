COMPILE_COMMANDS = -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -G Ninja
ARCHIVE = A2_SimpleNavigator_v1.0_CPP-1

all: install run

install:
	mkdir -p build
	cd build && cmake $(COMPILE_COMMANDS) .. && ninja

run:
	./build/console_app/main

console: s21_graph s21_graph_algorithms
	mkdir -p build
	cd build && cmake $(COMPILE_COMMANDS) .. && ninja console_app

s21_graph:
	mkdir -p build
	cd build && cmake $(COMPILE_COMMANDS) .. && ninja s21_graph

s21_graph_tests: s21_graph
	cd build && ninja graph_testing

s21_graph_algorithms:
	mkdir -p build
	cd build && cmake $(COMPILE_COMMANDS) .. && ninja s21_graph_algorithms

s21_graph_algorithms_tests: s21_graph_algorithms
	cd build && ninja graph_algo_testing

tests: s21_graph s21_graph_algorithms
	cd build && ninja graph_algo_testing && ninja graph_testing

uninstall:
	rm -rf build .cache

gcov_report: s21_graph s21_graph_algorithms 
	mkdir -p build
	mkdir -p report
	cd build && cmake $(COMPILE_COMMANDS) -DCOVERAGE=ON .. && ninja graph_algo_testing && ninja graph_testing && ninja coverage
	open report/report.html

dist:
	tar -cf $(ARCHIVE).tar console_app graph graph_algorithms lib tests CMakeLists.txt Makefile Doxyfile

dvi: 
	rm -rf docs
	mkdir docs
	doxygen Doxyfile
	mv html latex docs
	open docs/html/index.html

style:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -style=google -i ./source/graph/*.cpp ./source/graph/*.hpp ./source/graph_algorithms/*.cpp ./source/graph_algorithms/*.hpp ./console_app/*.cpp ./console_app/*.hpp
	clang-format -style=google -n ./source/graph/*.cpp ./source/graph/*.hpp ./source/graph_algorithms/*.cpp ./source/graph_algorithms/*.hpp ./console_app/*.cpp ./console_app/*.hpp
	rm .clang-format

clean: 
	rm -rf report/ docs/ build/ .cache/
