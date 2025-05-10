#include "s21_graph.hpp"

#include <ostream>

using namespace s21;

void Graph::ResizeMatrix(int size) {
  adjMatrix_.resize(size);
  for (auto &row : adjMatrix_) {
    row.resize(size, 0);
  }
}

bool Graph::LoadGraphFromFile(string filename) {
  ifstream f(filename);
  if (!f.is_open()) {
    cerr << "Error openning file: " << filename << endl;
    return false;
  }

  string line;
  getline(f, line);
  stringstream ss(line);

  if (!(ss >> nodesCount_) || nodesCount_ < 1) {
    cerr << "Error reading nodes count from file: " << filename << endl;
    return false;
  }

  ResizeMatrix(nodesCount_);
  int tmp;
  for (int i = 0; i < nodesCount_; ++i) {
    for (int j = 0; j < nodesCount_; ++j) {
      if (f >> tmp) {
        if (tmp < 0) {
          cerr << "Error! The matrix contains negative numbers!" << endl;
          return false;
        }
        if (tmp != 0 && tmp != 1) {
          isWeighted = true;
        }
        adjMatrix_[i][j] = tmp;
      } else {
        cerr << "Error reading matrix element at position (" << i << ", " << j
             << ")" << endl;
        return false;
      }
    }
  }
  return true;
}

bool Graph::IsDirected() {
  for (int i = 0; i < nodesCount_; ++i) {
    for (int j = 0; j < nodesCount_; ++j) {
      if (adjMatrix_[i][j] != adjMatrix_[j][i]) {
        return true;
      }
    }
  }
  return false;
}

bool Graph::ExportGraphToDot(string filename) {
  ofstream outfile(filename);
  if (!outfile.is_open()) {
    cerr << "Error opening file: " << filename << endl;
    return false;
  }

  bool directed = IsDirected();
  string nodesLinker = directed ? " -> " : " -- ";
  outfile << (directed ? "digraph" : "graph") << " graphname {" << endl;

  for (int i = 0; i < nodesCount_; ++i) {
    outfile << "  " << char(i + 97) << ";\n";
  }

  auto addEdge = [&](int i, int j) {
    if (adjMatrix_[i][j] != 0) {
      outfile << "  " << char(i + 97) << nodesLinker << char(j + 97);
      if (isWeighted) {
        outfile << " [weight=" << adjMatrix_[i][j] << "]";
      }
      outfile << ";\n";
    }
  };

  for (int i = 0; i < nodesCount_; ++i) {
    for (int j = directed ? 0 : i + 1; j < nodesCount_; ++j) {
      addEdge(i, j);
      if (!directed) {
        addEdge(j, i);
      }
    }
  }

  outfile << "}";
  outfile.close();
  return true;
}

const vector<vector<int>> &Graph::GetAdjMatrix() const noexcept {
  return adjMatrix_;
}
