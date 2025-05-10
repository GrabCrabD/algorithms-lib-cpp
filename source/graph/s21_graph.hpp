#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace s21 {

class Graph {
 public:
  Graph() : nodesCount_(0) {}

  Graph(const Graph &other)
      : nodesCount_(other.nodesCount_),
        adjMatrix_(other.adjMatrix_),
        isWeighted(other.isWeighted) {}

  Graph &operator=(const Graph &other) {
    if (this != &other) {
      nodesCount_ = other.nodesCount_;
      adjMatrix_ = other.adjMatrix_;
      isWeighted = other.isWeighted;
    }
    return *this;
  }

  ~Graph() {}

  /**
   * @brief Loading a graph from a file in the adjacency matrix format
   * @return true - ok, false - loading error
   */
  bool LoadGraphFromFile(string filename);

  /**
   * @brief Exporting a graph to a dot file
   * @return true - ok, false - exporting error
   */
  bool ExportGraphToDot(string filename);

  /**
   * @brief Getting graph's adjencency matrix
   * @return graph's adjencency matrix
   */
  const vector<vector<int>> &GetAdjMatrix() const noexcept;

  /**
   * @brief Checking graph is directed or not
   * @return true - directed, false - not directed
   */
  bool IsDirected();

 private:
  int nodesCount_;                 // adjacency matrix size
  vector<vector<int>> adjMatrix_;  // adjacency matrix
  bool isWeighted = false;         // weighted = true, unweighted = false

  void ResizeMatrix(int size);  // Resizing the matrix by size
};
}  // namespace s21
