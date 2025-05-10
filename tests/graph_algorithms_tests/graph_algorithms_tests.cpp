#include "s21_graph.hpp"
#include <gtest/gtest.h>

#include <s21_graph_algorithms.hpp>
#include <string>

using namespace s21;

TEST(S21GraphAlgorithmsTest, BreadthFirstSearch) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_FW.txt";
  g.LoadGraphFromFile(filename);

  GraphAlgorithms ga;
  vector<int> trueRes = {0, 1, 2, 3, 4};
  try {
    vector<int> res = ga.BreadthFirstSearch(g, 0);
    EXPECT_EQ(trueRes, res);
  } catch (string &ex) {
    FAIL();
  };
}

TEST(S21GraphAlgorithmsTest, BreadthFirstSearchEx) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_FW.txt";
  g.LoadGraphFromFile(filename);

  GraphAlgorithms ga;
  try {
    vector<int> res = ga.BreadthFirstSearch(g, 10);
  } catch (string &ex) {
    SUCCEED();
  };
}

TEST(S21GraphAlgorithmsTest, DepthFirstSearch) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_FW.txt";
  g.LoadGraphFromFile(filename);

  GraphAlgorithms ga;
  vector<int> trueRes = {0, 1, 2, 3, 4};
  try {
    vector<int> res = ga.DepthFirstSearch(g, 0);
    EXPECT_EQ(trueRes, res);
  } catch (string &ex) {
    FAIL();
  }
}

TEST(S21GraphAlgorithmsTest, DepthFirstSearchEx) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_FW.txt";
  g.LoadGraphFromFile(filename);

  GraphAlgorithms ga;
  try {
    vector<int> res = ga.DepthFirstSearch(g, 10);
  } catch (string &ex) {
    SUCCEED();
  };
}

TEST(S21GraphAlgorithmsTest, Dijkstra) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_FW.txt";
  g.LoadGraphFromFile(filename);

  GraphAlgorithms ga;
  try {
    int res = ga.GetShortestPathBetweenVertices(g, 0, 4);
    EXPECT_EQ(res, 8);
  } catch (string &ex) {
    FAIL();
  }
}

TEST(S21GraphAlgorithmsTest, FloydWarshall) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_FW.txt";
  g.LoadGraphFromFile(filename);

  GraphAlgorithms ga;
  vector<vector<int>> trueRes = {{0, 3, 9, 7, 8},
                                 {11, 0, 6, 4, 5},
                                 {5, 4, 0, 2, 4},
                                 {7, 2, 2, 0, 2},
                                 {10, 1, 5, 3, 0}};
  vector<vector<int>> res = ga.GetShortestPathsBetweenAllVertices(g);
  EXPECT_EQ(res, trueRes);
}

TEST(S21GraphAlgorithmsTest, AlgoPryme) {
  Graph g;
  string filename = "../../../dot_examples/adj_prim.txt";
  g.LoadGraphFromFile(filename);

  GraphAlgorithms ga;
  vector<vector<int>> trueRes = {{0, 9, 0, 0, 0},
                                 {9, 0, 0, 19, 0},
                                 {0, 0, 0, 51, 0},
                                 {0, 19, 51, 0, 31},
                                 {0, 0, 0, 31, 0}};
  try {
    vector<vector<int>> res = ga.GetLeastSpanningTree(g);
    EXPECT_EQ(trueRes, res);
  } catch (string &ex) {
    FAIL();
  }
}

TEST(S21GraphAlgorithmsTest, AlgoPrymeDirectGraphEx) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_FW.txt";
  g.LoadGraphFromFile(filename);

  GraphAlgorithms ga;
  try {
    vector<vector<int>> res = ga.GetLeastSpanningTree(g);
  } catch (string &ex) {
    SUCCEED();
  }
}

TEST(S21GraphAlgorithmsTest, AntAlgoMain) {
  Graph gCircle;
  GraphAlgorithms ga;

  string circleGraph = "../../../dot_examples/adj_matrix_ant_algo.txt";

  gCircle.LoadGraphFromFile(circleGraph);
  TsmResult tsmCircle = ga.SolveTravelingSalesmanProblem(gCircle);
  EXPECT_EQ(tsmCircle.distance, 40.0);
}

TEST(S21GraphAlgorithmsTest, AntAlgoDiff) {
  Graph gDiff;
  GraphAlgorithms ga;

  string difficultGraph = "../../../dot_examples/adj_matrix.txt";

  gDiff.LoadGraphFromFile(difficultGraph);
  TsmResult tsmDiff = ga.SolveTravelingSalesmanProblem(gDiff);
  if (tsmDiff.distance >= 120 && tsmDiff.distance <= 125) {
    SUCCEED();
  }
}

TEST(S21GraphAlgorithmsTest, GreedyAlgo) {
  Graph gCircle;
  GraphAlgorithms ga;

  string circleGraph = "../../../dot_examples/adj_matrix_ant_algo.txt";

  gCircle.LoadGraphFromFile(circleGraph);
  TsmResult tsmCircle = ga.TsmGreedyAlgo(gCircle);
  EXPECT_EQ(tsmCircle.distance, 40.0);
}

TEST(S21GraphAlgorithmsTest, GeneticAlgo) {
  Graph gCircle;
  GraphAlgorithms ga;

  string circleGraph = "../../../dot_examples/adj_matrix_ant_algo.txt";

  gCircle.LoadGraphFromFile(circleGraph);
  TsmResult tsmCircle = ga.GeneticAlgo(gCircle);
  EXPECT_EQ(tsmCircle.distance, 40.0);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
