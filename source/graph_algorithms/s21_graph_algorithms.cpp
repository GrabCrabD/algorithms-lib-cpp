#include "s21_graph_algorithms.hpp"

#include <sys/signal.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <ostream>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "s21_graph.hpp"
#include "s21_queue.hpp"

// genetic algorithm params:
constexpr int maxAttempts = 1000;  // max attempts of creating unique route
constexpr int generations = 100;   // count of max generations
constexpr int tournamentSize = 3;  // Tournament size

namespace s21 {

vector<int> GraphAlgorithms::BreadthFirstSearch(Graph &graph, int startVertex) {
  const vector<vector<int>> &matrix = graph.GetAdjMatrix();
  if (startVertex >= int(matrix.size()) || startVertex < 0) {
    throw string{"Error! Start vertex invalid!\n"};
  }

  Queue<int> que;
  vector<bool> visited(matrix.size(), false);
  vector<int> result;
  que.Push(startVertex);
  result.push_back(startVertex);
  visited[startVertex] = true;

  while (!que.Empty()) {
    int currentVertex = que.Front();
    que.Pop();
    for (size_t i = 0; i != matrix.size(); ++i) {
      if (matrix[currentVertex][i] != 0 && !visited[i]) {
        que.Push(i);
        visited[i] = true;
        result.push_back(i);
      }
    }
  }
  return result;
}

vector<int> GraphAlgorithms::DepthFirstSearch(Graph &graph, int startVertex) {
  const vector<vector<int>> &matrix = graph.GetAdjMatrix();
  if (startVertex >= int(matrix.size()) || startVertex < 0) {
    throw string{"Error! Start vertex invalid!\n"};
  }

  Stack<int> stack;
  vector<int> result;
  vector<bool> visited(matrix.size(), false);

  stack.Push(startVertex);
  result.push_back(startVertex);
  visited[startVertex] = true;

  while (!stack.Empty()) {
    int currentVertex = stack.Top();
    stack.Pop();
    for (size_t i = 0; i != matrix.size(); ++i) {
      if (matrix[currentVertex][i] != 0 && !visited[i]) {
        stack.Push(i);
        visited[i] = true;
        result.push_back(i);
      }
    }
  }
  return result;
}

int GraphAlgorithms::GetShortestPathBetweenVertices(Graph &graph, int vertex1,
                                                    int vertex2) {
  const vector<vector<int>> &matrix = graph.GetAdjMatrix();
  int matrixSize = matrix.size();
  if (vertex2 >= matrixSize) {
    cout << "vertex2 = " << vertex2 << " matrix size = " << matrixSize << endl;
    throw string{"Error! Vertex2 invalid!\n"};
  } else if (vertex1 < 0 || vertex1 >= matrixSize) {
    throw string{"Error! Vertex1 invalid!\n"};
  }

  Queue<int> que;
  vector<int> distance(matrixSize, numeric_limits<int>::max());
  vector<bool> visited(matrixSize, false);

  distance[vertex1] = 0;
  que.Push(vertex1);

  while (!que.Empty()) {
    int currentVertex = que.Front();
    que.Pop();

    if (visited[currentVertex]) {
      continue;
    }
    visited[currentVertex] = true;

    for (int i = 0; i != matrixSize; ++i) {
      if (matrix[currentVertex][i] != 0 && !visited[i]) {
        int newDist = distance[currentVertex] + matrix[currentVertex][i];
        if (newDist < distance[i]) {
          distance[i] = newDist;
          que.Push(i);
        }
      }
    }
  }

  if (distance[vertex2] == numeric_limits<int>::max()) {
    throw string{"No paths from vertex1 to vertex2!\n"};
  }

  return distance[vertex2];
}

vector<vector<int>> GraphAlgorithms::GetShortestPathsBetweenAllVertices(
    Graph &graph) {
  const vector<vector<int>> &matrix = graph.GetAdjMatrix();
  int matrix_size = matrix.size();

  vector<vector<int>> distances(matrix_size, vector<int>(matrix_size));

  for (int i = 0; i < matrix_size; ++i) {
    for (int j = 0; j != matrix_size; ++j) {
      // filling every 0 in main matrix with inf but values
      distances[i][j] =
          matrix[i][j] == 0 ? numeric_limits<int>::max() : matrix[i][j];
      // filling the main diagonal with 0
      distances[i][j] = i == j ? 0 : distances[i][j];
    }
  }

  for (int v = 0; v < matrix_size; ++v) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < matrix_size; ++j) {
        if (distances[i][v] != numeric_limits<int>::max() &&
            distances[v][j] != numeric_limits<int>::max() &&
            distances[i][j] > distances[i][v] + distances[v][j]) {
          distances[i][j] = distances[i][v] + distances[v][j];
        }
      }
    }
  }

  return distances;
}

vector<vector<int>> GraphAlgorithms::GetLeastSpanningTree(Graph &graph) {
  if (graph.IsDirected()) {
    throw string{"Error! Graph is directed!\n"};
  }

  const vector<vector<int>> &matrix = graph.GetAdjMatrix();
  size_t matrixSize = matrix.size();
  vector<vector<int>> spanTreeMatrix(matrixSize, vector<int>(matrixSize, 0));
  vector<bool> visited(matrixSize, false);
  vector<int> parent(matrixSize, -1);
  vector<int> minEdge(matrixSize, numeric_limits<int>::max());

  minEdge[0] = 0;

  for (size_t i = 0; i != matrixSize; ++i) {
    int minKey = numeric_limits<int>::max();
    int minIdx = -1;

    // Serach for vertex with min weight
    for (size_t j = 0; j != matrixSize; ++j) {
      if (!visited[j] && minEdge[j] < minKey) {
        minKey = minEdge[j];
        minIdx = j;
      }
    }

    visited[minIdx] = true;

    // Update minimum edges for neighboring vertices
    for (size_t k = 0; k != matrixSize; ++k) {
      if (matrix[minIdx][k] != 0 && !visited[k] &&
          matrix[minIdx][k] < minEdge[k]) {
        parent[k] = minIdx;
        minEdge[k] = matrix[minIdx][k];
      }
    }
  }

  // Filling matrix MST
  for (size_t i = 1; i != matrixSize; ++i) {
    if (parent[i] != -1) {
      spanTreeMatrix[parent[i]][i] = matrix[i][parent[i]];
      spanTreeMatrix[i][parent[i]] = matrix[i][parent[i]];  // undirected graphs
    }
  }

  return spanTreeMatrix;
}

TsmResult GraphAlgorithms::SolveTravelingSalesmanProblem(Graph &graph) {
  const vector<vector<int>> &matrix = graph.GetAdjMatrix();
  const int m = matrix.size();  // vertex's count

  TsmResult tsm;
  // size = m+1 to add path from last node to start
  tsm.vertices.resize(m + 1, 0);
  tsm.distance = 0.0;

  int exitIndex = 0;  // index that count kExitCondition last iteration

  // data of the minimal distances of the last kExitCondition iteration
  vector<double> exitCondition(kExitCondition, 0.0);

  // global pheromone matrix data
  vector<vector<double>> globalPheromoneData(m, vector(m, kPheromoneStart));

  // local pheromone matrix data
  vector<vector<double>> localPheromoneData(m, vector(m, 0.0));

  for (int current = 0; current != iterationCount; ++current) {
    // path distance that walked each ant for one iteration
    vector<double> antDistance(m, 0.0);

    // path's nodes walked by each ant for one iteration
    vector<vector<int>> antPath(m, vector(m + 1, 0));

    for (int i = 0; i != m; ++i) {
      oneIteration(i, globalPheromoneData, matrix, m, alpha, beta, antDistance,
                   antPath);
    }

    for (int i = 0; i != m; ++i) {
      localPheromoneUpdate(localPheromoneData, antDistance[i], antPath[i], Q,
                           m);
    }

    globalPheromoneUpdate(globalPheromoneData, localPheromoneData, m,
                          kPheromoneVaporizing);

    tsmUpdate(antDistance, antPath, tsm);

    if (!exitConditionCheck(exitIndex, kExitCondition, exitCondition, tsm)) {
      break;
    }
  }

  return tsm;
}

bool GraphAlgorithms::exitConditionCheck(int &exitIndex,
                                         const int kExitCondition,
                                         vector<double> &exitCondition,
                                         TsmResult &tsm) {
  // checking for filling exitCondition vector
  if (exitIndex != kExitCondition) {
    exitCondition[exitIndex] = tsm.distance;
    exitIndex++;
  } else {
    exitIndex = 0;
    exitCondition[exitIndex] = tsm.distance;
  }

  double exitSum = 0.0;  // checking for repeatable minimal distance
  for (int i = 0; i != kExitCondition; ++i) {
    exitSum += exitCondition[i];
  }
  if (exitCondition[0] == exitSum / kExitCondition) {
    return false;
  }
  return true;
}

void GraphAlgorithms::localPheromoneUpdate(
    vector<vector<double>> &localPheromoneData, double &antDistance,
    vector<int> &antPath, const double &Q, const int &m) {
  for (int i = 0; i != m - 1; ++i) {
    // the sum of local pheromone updates for the current vertex
    localPheromoneData[antPath[i]][antPath[i + 1]] += Q / antDistance;
  }
}

void GraphAlgorithms::globalPheromoneUpdate(
    vector<vector<double>> &globalPheromoneData,
    vector<vector<double>> &localPheromoneData, const int &m,
    const double &kPheromoneVaporizing) {
  for (int i = 0; i != m; ++i) {
    for (int j = 0; j != m; ++j) {
      globalPheromoneData[i][j] =
          (globalPheromoneData[i][j] * kPheromoneVaporizing) +
          localPheromoneData[i][j];
    }
  }
}

void GraphAlgorithms::tsmUpdate(vector<double> &antDistance,
                                vector<vector<int>> antPath, TsmResult &tsm) {
  auto it = std::min_element(antDistance.begin(), antDistance.end());
  double min = *it;
  int iMin = std::distance(antDistance.begin(), it);

  if (tsm.distance == 0.0 || tsm.distance >= min) {
    tsm.distance = min;
    tsm.vertices = antPath[iMin];
  }
}

void GraphAlgorithms::oneIteration(int antNum,
                                   vector<vector<double>> &globalPheromoneData,
                                   const vector<vector<int>> &matrix,
                                   const int &m, const int &alpha,
                                   const int &beta, vector<double> &antDistance,
                                   vector<vector<int>> &antPath) {
  int current = antNum;
  int currentTmp = 0;

  vector<bool> visited(m, false);  // data of visited nodes by ants

  for (int i = 0; i != m - 1; ++i) {  // ant's steps
    antPath[antNum][i] = current;
    visited[current] = true;
    currentTmp = current;  // current node's index
    vector<double> prob = probabilityCalculation(
        current, globalPheromoneData, matrix, m, alpha, beta, visited);
    current = stepCalculation(prob);  // next node's index

    antDistance[antNum] += matrix[currentTmp][current];
  }

  antPath[antNum][m - 1] = current;
  antPath[antNum][m] = antNum;

  antDistance[antNum] += matrix[current][antNum];
}

vector<double> GraphAlgorithms::probabilityCalculation(
    int current, vector<vector<double>> &globalPheromoneData,
    const vector<vector<int>> &matrix, const int &m, const int &alpha,
    const int &beta, vector<bool> &visited) {
  vector<double> antWish(m, 0);  // vector of ant's wishes in one step

  for (int i = 0; i != m; ++i) {
    if (!visited[i]) {
      antWish[i] = pow(globalPheromoneData[current][i], alpha) *
                   pow(100 / matrix[current][i], beta);
    }
  }

  double sumWishes = 0;
  for (auto &antWish : antWish) {
    sumWishes += antWish;
  }

  vector<double> probability(m, 0);
  for (int i = 0; i != m; ++i) {
    probability[i] = (antWish[i] / sumWishes) * 100;
  }
  return probability;
}

int GraphAlgorithms::stepCalculation(vector<double> stepProbabilities) {
  double rand = random_index();
  double sum = 0;
  int nextStep = 0;
  for (int i = 0; rand >= sum; ++i) {
    sum += stepProbabilities[i];
    nextStep = i;
  }
  return nextStep;
}

double GraphAlgorithms::random_index() {
  std::random_device rd;
  std::default_random_engine engine(rd());
  std::uniform_real_distribution<double> dist(0.0, 100.0);
  return dist(engine);
}

TsmResult GraphAlgorithms::TsmGreedyAlgo(Graph &graph) {
  TsmResult tsm;
  const vector<vector<int>> &matrix = graph.GetAdjMatrix();
  size_t matrixSize = matrix.size();
  vector<bool> visited(matrixSize, false);

  srand(time(0));

  // init first vertex
  int currentVertex = rand() % matrixSize;
  visited[currentVertex] = true;
  tsm.vertices.push_back(currentVertex);

  for (size_t step = 1; step < matrixSize; ++step) {
    int nextVertex = -1;
    int minEdge = numeric_limits<int>::max();

    for (size_t i = 0; i < matrixSize; ++i) {
      if (!visited[i] && matrix[currentVertex][i] &&
          matrix[currentVertex][i] < minEdge) {
        minEdge = matrix[currentVertex][i];
        nextVertex = i;
      }
    }

    if (nextVertex == -1) {
      throw logic_error("Error! Graph is not fully connected.\n");
    }

    // refresh current vertex and distance
    currentVertex = nextVertex;
    tsm.vertices.push_back(currentVertex);
    tsm.distance += minEdge;
    visited[currentVertex] = true;
  }

  // add last distance between start and last vertex
  int returnDistance = matrix[currentVertex][tsm.vertices[0]];
  if (returnDistance != 0) {
    tsm.distance += returnDistance;
    tsm.vertices.push_back(tsm.vertices[0]);
  } else {
    throw logic_error("No return path to the starting vertex.\n");
  }

  return tsm;
}

TsmResult GraphAlgorithms::GeneticAlgo(Graph &graph) {
  TsmResult tsm;
  const vector<vector<int>> &matrix = graph.GetAdjMatrix();
  size_t populationSize = matrix.size();
  size_t numCities = matrix.size();

  // create initial population
  vector<vector<int>> population =
      CreateInitialPopulation(populationSize, numCities, matrix);

  int noChangeCount = 0;
  TsmResult bestResult;

  for (int gen = 0; gen != generations; ++gen) {
    // calculate fitness for current population
    vector<int> fitness = CalculateFitness(population, matrix);

    // buffer population for parents and descendants
    vector<vector<int>> bufferPopulation = population;

    // tournament choosing pair of parents
    while (bufferPopulation.size() < populationSize * 2) {
      pair<vector<int>, vector<int>> parents =
          TournamentSelection(population, fitness);

      pair<vector<int>, vector<int>> descendants =
          Crossbreeding(parents.first, parents.second);

      Mutation(descendants.first);
      Mutation(descendants.second);

      bufferPopulation.push_back(descendants.first);
      bufferPopulation.push_back(descendants.second);
    }

    // fitness for bufferPopulation (parents and descendants)
    vector<int> bufferFitness = CalculateFitness(bufferPopulation, matrix);

    vector<TsmResult> bufPopulationWithFitness(bufferPopulation.size());
    for (size_t i = 0; i != bufferPopulation.size(); ++i) {
      bufPopulationWithFitness[i].vertices = bufferPopulation[i];
      bufPopulationWithFitness[i].distance = bufferFitness[i];
    }

    sort(bufPopulationWithFitness.begin(), bufPopulationWithFitness.end(),
         [](const TsmResult &first, const TsmResult &second) {
           return first.distance < second.distance;
         });

    vector<vector<int>> newPopulation;
    for (size_t i = 0; i != bufPopulationWithFitness.size() &&
                       newPopulation.size() != populationSize;
         ++i) {
      newPopulation.push_back(bufPopulationWithFitness[i].vertices);
    }

    population = newPopulation;

    // tsm.vertices = bufPopulationWithFitness[0].vertices;
    // tsm.distance = bufPopulationWithFitness[0].distance;

    if (gen == 0 ||
        bufPopulationWithFitness[0].distance < bestResult.distance) {
      bestResult = bufPopulationWithFitness[0];
      noChangeCount = 0;
    } else {
      noChangeCount++;
    }

    if (noChangeCount > kExitCondition) {
      break;
    }
  }

  if (!bestResult.vertices.empty()) {
    bestResult.vertices.push_back(bestResult.vertices[0]);
  }

  tsm = bestResult;

  return tsm;
}

vector<vector<int>> GraphAlgorithms::CreateInitialPopulation(
    size_t populationSize, size_t numCities,
    const vector<vector<int>> &matrix) {
  vector<vector<int>> population;
  set<vector<int>> uniqueRoutes;

  // init base route
  vector<int> baseRoute(numCities);
  for (size_t i = 0; i != numCities; ++i) {
    baseRoute[i] = i;
  }

  random_device rd;
  mt19937 g(rd());
  int attempts = 0;

  while (population.size() < populationSize && attempts < maxAttempts) {
    vector<int> newRoute = baseRoute;
    shuffle(newRoute.begin(), newRoute.end(), g);

    // checking the route for uniqueness and validity
    if (uniqueRoutes.find(newRoute) == uniqueRoutes.end() &&
        IsValidRoute(newRoute, matrix)) {
      uniqueRoutes.insert(newRoute);
      population.push_back(newRoute);
      attempts = 0;
    } else {
      attempts++;
    }
  }

  if (population.size() < populationSize) {
    cerr << "Warning: Could not generate the required number of unique "
            "routes.\n";
    cerr << "Try to increase maxAttrmpts value in function: " << (__FUNCTION__)
         << " or load new graph.\n";
  }

  return population;
}

bool GraphAlgorithms::IsValidRoute(const vector<int> &route,
                                   const vector<vector<int>> &matrix) {
  // Check if there's a valid path between each consecutive city in the route
  for (size_t i = 0; i < route.size() - 1; ++i) {
    if (matrix[route[i]][route[i + 1]] == 0) {  // No path between cities
      return false;
    }
  }

  // Also check the path returning to the starting city
  if (matrix[route.back()][route[0]] == 0) {
    return false;
  }

  return true;
}

pair<vector<int>, vector<int>> GraphAlgorithms::Crossbreeding(
    const vector<int> &parentX, const vector<int> &parentY) {
  random_device rd;
  mt19937 g(rd());
  uniform_int_distribution<> distrib(1, parentX.size() - 1);
  int alpha = distrib(g);  // init a crossover point

  vector<int> descendantX(parentX.size()), descendantY(parentY.size());
  // copying the first part of each parent to their  descendants
  for (int i = 0; i < alpha; ++i) {
    descendantX[i] = parentX[i];
    descendantY[i] = parentY[i];
  }

  // filling out the other part of each descendant
  auto FillRemaining = [](const vector<int> &parent, vector<int> &descendant,
                          int alpha) {
    size_t j = alpha;
    for (size_t i = 0; i != parent.size(); ++i) {
      if (find(descendant.begin(), descendant.end(), parent[i]) ==
          descendant.end()) {
        descendant[j++] = parent[i];
      }
    }
  };

  FillRemaining(parentX, descendantY, alpha);
  FillRemaining(parentY, descendantX, alpha);

  return make_pair(descendantX, descendantY);
}

void GraphAlgorithms::Mutation(vector<int> &descendant) {
  random_device rd;
  mt19937 g(rd());

  // mutation probability 2 / chromosome_length
  double mutationProbability = 0.5;
  uniform_real_distribution<> probDistrib(0.0, 1.0);

  if (probDistrib(g) < mutationProbability) {
    uniform_int_distribution<> startDistrib(0, descendant.size() - 2);
    int startGen = startDistrib(g);

    uniform_int_distribution<> endDistrib(startGen + 1, descendant.size() - 1);
    int endGen = endDistrib(g);
    reverse(descendant.begin() + startGen, descendant.begin() + endGen + 1);
  }
}

vector<int> GraphAlgorithms::CalculateFitness(
    const vector<vector<int>> &population, const vector<vector<int>> &matrix) {
  vector<int> fitness(population.size());

  for (size_t i = 0; i != population.size(); ++i) {
    for (size_t j = 0; j != population[i].size() - 1; ++j) {
      int from = population[i][j];
      int to = population[i][j + 1];
      fitness[i] += matrix[from][to];
    }

    // plus distance between start and end cities
    fitness[i] += matrix[population[i].back()][population[i].front()];
  }

  return fitness;
}

pair<vector<int>, vector<int>> GraphAlgorithms::TournamentSelection(
    const vector<vector<int>> &population, const vector<int> &fitness) {
  random_device rd;
  mt19937 g(rd());
  uniform_int_distribution<> distrib(0, population.size() - 1);

  // choosing first parent
  vector<int> tournament1(tournamentSize);
  for (size_t i = 0; i != tournamentSize; ++i) {
    tournament1[i] = distrib(g);
  }

  int bestIndex1 = tournament1[0];
  for (size_t i = 1; i != tournamentSize; ++i) {
    if (fitness[tournament1[i]] < fitness[bestIndex1]) {
      bestIndex1 = tournament1[i];
    }
  }
  vector<int> parentX = population[bestIndex1];

  // choosing second unique parent
  vector<int> tournament2(tournamentSize);
  int bestIndex2;
  do {
    for (size_t i = 0; i != tournamentSize; ++i) {
      tournament2[i] = distrib(g);
    }

    bestIndex2 = tournament2[0];
    for (size_t i = 1; i != tournamentSize; ++i) {
      if (fitness[tournament2[i]] < fitness[bestIndex2]) {
        bestIndex2 = tournament2[i];
      }
    }
  } while (bestIndex2 == bestIndex1);  // checking for unique

  vector<int> parentY = population[bestIndex2];

  return make_pair(parentX, parentY);
}

}  // namespace s21
