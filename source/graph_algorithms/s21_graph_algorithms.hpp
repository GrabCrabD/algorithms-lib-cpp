#pragma once

#include <cxxabi.h>
#include <sys/signal.h>

#include <cstddef>
#include <cstdlib>
#include <s21_graph.hpp>
#include <s21_queue.hpp>
#include <s21_stack.hpp>
#include <vector>

namespace s21 {

struct TsmResult {
  vector<int> vertices;   // an array with the route
  double distance = 0.0;  // the length of this route
};

// variables for AntAlgo

constexpr int alpha = 1;                      // ratio pheromone
constexpr int beta = 1;                       // ratio distance
constexpr int iterationCount = 100;           // iteration count
constexpr double kPheromoneStart = 1;         // strart pheromone count
constexpr double kPheromoneVaporizing = 0.8;  // pheromone vaporizing ratio
constexpr double Q = 101;                     // ranged constant
constexpr int kExitCondition =
    10;  // count of last iteration for exit condition

class GraphAlgorithms {
 public:
  /**
   * @brief a non-recursive depth-first search in the graph from a given vertex
   * @return an array that contains the traversed vertices in the order they
   * were traversed
   */
  vector<int> DepthFirstSearch(Graph &graph, int start_vertex);

  /**
   * @brief breadth-first search in the graph from a given vertex
   * @return an array that contains the traversed vertices in the order they
   * were traversed
   */
  vector<int> BreadthFirstSearch(Graph &graph, int start_vertex);

  /**
   * @brief searching for the shortest path between two vertices in a graph
   * using Dijkstra's algorithm
   * @param graph graph object
   * @param vertex1 start vertex
   * @param vertex2 end vertex
   * @return the smallest distance between vertices
   */
  int GetShortestPathBetweenVertices(Graph &graph, int vertex1, int vertex2);

  /**
   * @brief searching for the shortest paths between all pairs of vertices in a
   * graph using the Floyd-Warshall algorithm
   * @param graph graph object
   * @return the matrix of the shortest paths between all vertices of the graph
   */
  vector<vector<int>> GetShortestPathsBetweenAllVertices(Graph &graph);

  /**
   * @brief searching for the minimal spanning tree in a graph using Prim's
   * algorithm
   * @param graph graph object
   * @return the adjacency matrix for the minimal spanning tree
   */
  vector<vector<int>> GetLeastSpanningTree(Graph &graph);

  /**
   * @brief solving the traveling salesman's problem using the ant colony
   * algorithm. Find the shortest path that goes through all vertices of the
   * graph at least once, followed by a return to the original vertex
   * @param graph graph object
   * @return result data of solving the traveling salesman's problem
   */
  TsmResult SolveTravelingSalesmanProblem(Graph &graph);

  /**
   * @brief solving the traveling salesman's problem using the greedy
   * algorithm. Find the shortest path that goes through all vertices of the
   * graph at least once, followed by a return to the original vertex
   * @return the `TsmResult` structure
   */
  TsmResult TsmGreedyAlgo(Graph &graph);

  /**
   * @brief solving the traveling salesman's problem using the genetic algorithm
   * @param graph graph object
   */
  TsmResult GeneticAlgo(Graph &graph);

 private:
  /**
  @brief solving one iteration of the traveling salesman's problem
  @param antNum ant number
  @param globalPheromoneData global pheromone matrix
  @param matrix graph's adjency matrix
  @param m graph's nodes count
  @param alpha pheromone ratio
  @param beta distance ratio
  @param antDistance distance data for each ant in one iteration
  @param antPath path data for each ant in one iteration
   */
  void oneIteration(int antNum, vector<vector<double>> &globalPheromoneData,
                    const vector<vector<int>> &matrix, const int &m,
                    const int &alpha, const int &beta,
                    vector<double> &antDistance, vector<vector<int>> &antPath);

  /**
  @brief calculating probability for each ant for the next step in one iteration
  @param start current node for each ant
  @param globalPheromoneData global pheromone matrix
  @param matrix graph's adjency matrix
  @param m graph's nodes count
  @param alpha pheromone ratio
  @param beta distance ratio
  @param visited data of visited nodes
  @return data of probabilities to step for each nodes
   */
  vector<double> probabilityCalculation(
      int current, vector<vector<double>> &globalPheromoneData,
      const vector<vector<int>> &matrix, const int &m, const int &alpha,
      const int &beta, vector<bool> &visited);

  /**
  @brief calculating local pheromone updating in one iteration
  @param localPheromoneData local pheromone matrix
  @param antDistance distance data for each ant in one iteration
  @param antPath path data for each ant in one iteration
  @param Q local pheromone updating ratio
  @param m graph's nodes count
   */
  void localPheromoneUpdate(vector<vector<double>> &localPheromoneData,
                            double &antDistance, vector<int> &antPath,
                            const double &Q, const int &m);

  /**
  @brief calculating global pheromone updating in one iteration
  @param globalPheromoneData global pheromone matrix
  @param localPheromoneData local pheromone matrix
  @param m graph's nodes count
  @param kPheromoneVaporizing pheromone vaporizing ratio
  */
  void globalPheromoneUpdate(vector<vector<double>> &globalPheromoneData,
                             vector<vector<double>> &localPheromoneData,
                             const int &m, const double &kPheromoneVaporizing);

  /**
  @brief updating result data of solving the traveling salesman's problem
  @param antDistance distance data for each ant in one iteration
  @param antPath path data for each ant in one iteration
  @param tsm structure with result data of solving the traveling salesman's
  problem
  */
  void tsmUpdate(vector<double> &antDistance, vector<vector<int>> antPath,
                 TsmResult &tsm);

  /**
  @brief checking exit conditions from solving the traveling salesman's problem
  @param exitIndex index that count kExitCondition last iteration
  @param kExitCondition count of last iteration for exit condition
  @param exitCondition data of the minimal distances of the last kExitCondition
  iteration

  */
  bool exitConditionCheck(int &exitIndex, const int kExitCondition,
                          vector<double> &exitCondition, TsmResult &tsm);

  /**
  @brief calculation to which next node step
  @param stepProbabilities probabilities data for each node for the next
  step
  @return node's number for the next step
  */
  int stepCalculation(vector<double> stepProbabilities);

  double random_index();

  vector<vector<int>> CreateInitialPopulation(
      size_t populationSize, size_t numCities,
      const vector<vector<int>> &matrix);
  pair<vector<int>, vector<int>> Crossbreeding(const vector<int> &parentX,
                                               const vector<int> &parentY);

  vector<int> CalculateFitness(const vector<vector<int>> &population,
                               const vector<vector<int>> &matrix);
  void Mutation(vector<int> &descendant);

  pair<vector<int>, vector<int>> TournamentSelection(
      const vector<vector<int>> &population, const vector<int> &fitness);

  bool IsValidRoute(const vector<int> &route,
                    const vector<vector<int>> &matrix);
};

}  // namespace s21
