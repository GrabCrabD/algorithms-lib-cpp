#include <gtest/gtest.h>

#include <s21_graph.hpp>

using namespace s21;

TEST(S21GraphTest, LoadGraphFromFileFailure1) {
  Graph g;
  string filename = "failure.txt";

  ASSERT_FALSE(g.LoadGraphFromFile(filename));
  ASSERT_FALSE(g.ExportGraphToDot("/?s.txt"));
}

TEST(S21GraphTest, LoadGraphFromFileFailure2) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_error_file.txt";

  ASSERT_FALSE(g.LoadGraphFromFile(filename));
}

TEST(S21GraphTest, LoadGraphFromFileFailure3) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_error.txt";

  ASSERT_FALSE(g.LoadGraphFromFile(filename));
}

TEST(S21GraphTest, LoadGraphFromFileSuccess) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix1.txt";

  ASSERT_TRUE(g.LoadGraphFromFile(filename));
}

TEST(S21GraphTest, LoadAndExportGraphFromFileWeighted) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_weighted.txt";

  ASSERT_TRUE(g.LoadGraphFromFile(filename));
  ASSERT_TRUE(g.ExportGraphToDot("wei_result.txt"));
}

TEST(S21GraphTest, LoadAndExportGraphDirected) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_directed.txt";

  ASSERT_TRUE(g.LoadGraphFromFile(filename));
  ASSERT_TRUE(g.ExportGraphToDot("dir_result.txt"));
}

TEST(S21GraphTest, LoadAndExportGraphDirectedWeigted) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_weighted.txt";

  ASSERT_TRUE(g.LoadGraphFromFile(filename));
  ASSERT_TRUE(g.ExportGraphToDot("dir_wei_result.txt"));
}

TEST(S21GraphTest, LoadAndExportGraphFromFile) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix1.txt";

  ASSERT_TRUE(g.LoadGraphFromFile(filename));
  ASSERT_TRUE(g.ExportGraphToDot("matrix1_result.txt"));
}

TEST(S21GraphTest, LoadAndExportGraphFromFileLoop) {
  Graph g;
  string filename = "../../../dot_examples/adj_matrix_loop.txt";

  ASSERT_TRUE(g.LoadGraphFromFile(filename));
  ASSERT_TRUE(g.ExportGraphToDot("loop_result.txt"));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
