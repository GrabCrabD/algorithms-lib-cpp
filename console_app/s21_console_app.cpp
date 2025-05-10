#include "s21_console_app.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

#include "s21_graph_algorithms.hpp"

namespace s21 {
bool ConsoleApp::SetFileNameFromConsole() {
  bool status = false;
  char *line = readline("Enter filename: ");

  if (!line) {
    throw "Error! Input line invalid\n";
  }
  if (*line) {
    add_history(line);
  }

  std::string filename = std::string(line);

  if (g_.LoadGraphFromFile(filename)) {
    std::cout << "Graph loaded successfully.\n";
    status = true;
  }
  free(line);
  return status;
}

void ConsoleApp::ConsoleBFS() {
  try {
    vector<int> result = ga_.BreadthFirstSearch(g_, vertex1_);
    for (size_t i = 0; i != result.size(); ++i) {
      cout << " " << result[i] << " ";
    }
  } catch (string &errorMessage) {
    cerr << errorMessage;
  }
}

void ConsoleApp::ConsoleDFS() {
  try {
    vector<int> result = ga_.DepthFirstSearch(g_, vertex1_);
    for (size_t i = 0; i != result.size(); ++i) {
      cout << " " << result[i] << " ";
    }
  } catch (string &errorMessage) {
    cerr << errorMessage;
  }
}

void ConsoleApp::ConsoleDijkstra() {
  int res;
  try {
    res = ga_.GetShortestPathBetweenVertices(g_, vertex1_, vertex2_);
    cout << "Shortest path between vertices = " << res << endl;
  } catch (string &errorMessage) {
    cerr << errorMessage;
  }
}

void ConsoleApp::ConsoleFloydWarshall() {
  vector<vector<int>> distancesMatrix =
      ga_.GetShortestPathsBetweenAllVertices(g_);

  for (size_t i = 0; i != distancesMatrix.size(); ++i) {
    for (size_t k = 0; k != distancesMatrix.size(); ++k) {
      cout << " " << distancesMatrix[i][k] << " ";
    }
    cout << endl;
  }
}

void ConsoleApp::ConsolePrim() {
  try {
    vector<vector<int>> distancesMatrix = ga_.GetLeastSpanningTree(g_);

    for (size_t i = 0; i != distancesMatrix.size(); ++i) {
      for (size_t k = 0; k != distancesMatrix.size(); ++k) {
        cout << " " << distancesMatrix[i][k] << " ";
      }
      cout << endl;
    }
  } catch (string &errorMessage) {
    cerr << errorMessage;
  }
}

void ConsoleApp::ConsoleAntAlgo() {
  TsmResult tsmR = ga_.SolveTravelingSalesmanProblem(g_);

  cout << "TsmResult:\n";
  cout << "vertices:";
  for (size_t i = 0; i != tsmR.vertices.size(); ++i) {
    cout << " " << tsmR.vertices[i] << " ";
  }
  cout << "\ndistance = " << tsmR.distance << endl;
}

void ConsoleApp::ComparisonTests() {
  TsmResult tsmR1, tsmR2, tsmR3;
  cout << '\n';

  {
    LogDuration logDur("Ant algo");
    for (int i = 0; i != N_; ++i) {
      try {
        tsmR1 = ga_.SolveTravelingSalesmanProblem(g_);
      } catch (const std::logic_error &e) {
        cerr << e.what();
        cerr << "Load new graph or choose other menu option\n";
        return;
      }
    }
  }

  for (size_t i = 0; i != tsmR1.vertices.size(); ++i) {
    cout << " " << tsmR1.vertices[i] << " ";
  }
  cout << "\ndistance = " << tsmR1.distance << '\n' << '\n';

  {
    LogDuration logDur("Greedy algo");
    for (int i = 0; i != N_; ++i) {
      try {
        tsmR2 = ga_.TsmGreedyAlgo(g_);
      } catch (const std::logic_error &e) {
        cerr << e.what();
        cerr << "Load new graph or choose other menu option\n";
        return;
      }
    }
  }

  for (size_t i = 0; i != tsmR2.vertices.size(); ++i) {
    cout << " " << tsmR2.vertices[i] << " ";
  }
  cout << "\ndistance = " << tsmR2.distance << '\n' << '\n';

  {
    LogDuration logDur("Genetic algo");
    for (int i = 0; i != N_; ++i) {
      try {
        tsmR3 = ga_.GeneticAlgo(g_);
      } catch (const std::logic_error &e) {
        cerr << e.what();
        cerr << "Load new graph or choose other menu option\n";
        return;
      }
    }
  }

  for (size_t i = 0; i != tsmR3.vertices.size(); ++i) {
    cout << " " << tsmR3.vertices[i] << " ";
  }
  cout << "\ndistance = " << tsmR3.distance << '\n';
}

void ConsoleApp::DisplayMenu() {
  cout << "\n==== Graph Manager ====\n";
  cout << "1. Load new graph from file\n";
  cout << "2. Traverse the graph in breadth and print the result\n";
  cout << "3. Traverse the graph in depth and print the result\n";
  cout << "4. Find the shortest path between any two vertices and print the "
          "result\n";
  cout << "5. Find the shortest paths between all pairs of vertices in the "
          "graph and print the result matrix\n";
  cout << "6. Search for the minimum spanning tree in the graph and print the "
          "resulting adjacency matrix\n";
  cout << "7. Solve the Salesman problem, with output of the resulting route "
          "and its length\n";
  cout << "8. Bonus. Comparison of methods for solving the traveling salesman "
          "problem\n";
  cout << "0. Exit\n";
  cout << "=======================\n";
  cout << "Enter your choice: ";
}

void ConsoleApp::SetVertices(int vertex1, int vertex2) {
  vertex1_ = vertex1;
  vertex2_ = vertex2;
}

void ConsoleApp::SetStartVertex(int startVertex) { vertex1_ = startVertex; }

void ConsoleApp::SetN(int N) { N_ = N; }
}  // namespace s21
