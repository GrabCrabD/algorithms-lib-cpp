#pragma once

#include <readline/history.h>
#include <readline/readline.h>

#include <s21_graph_algorithms.hpp>

namespace s21 {
class ConsoleApp {
 public:
  void DisplayMenu();
  bool SetFileNameFromConsole();
  void ConsoleDijkstra();
  void ConsoleBFS();
  void ConsoleDFS();
  void ConsoleFloydWarshall();
  void ConsolePrim();
  void ConsoleAntAlgo();  // salesman problem; part 4
  void ComparisonTests();

  void SetVertices(int vertex1, int vertex2);
  void SetStartVertex(int start_vertex);
  void SetN(int N);

  class LogDuration {
   public:
    LogDuration(string id) : id_(std::move(id)) {}

    ~LogDuration() {
      const auto end_time = chrono::steady_clock::now();
      const auto dur = end_time - start_time_;
      cout << id_ << ": ";
      cout << "operation time: "
           << chrono::duration_cast<chrono::milliseconds>(dur).count() << " ms"
           << endl;
    }

   private:
    const string id_;
    const chrono::steady_clock::time_point start_time_ =
        chrono::steady_clock::now();
  };

 private:
  int vertex1_, vertex2_;
  int N_;
  Graph g_;
  GraphAlgorithms ga_;
};
}  // namespace s21
