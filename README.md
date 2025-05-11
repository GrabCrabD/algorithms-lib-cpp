# Simple Navigator

## Build

  ```bash
  # Build libraries and tests
  make

  # Run tests
  make tests

  # Create doxygen docs
  make dvi

  # Create tar archive
  make dist

  # Check tests coverege
  make gcov_report

  # Clean
  make clean
  ```

## Key Features

### Core Components

- **Graph Representation**  
  Dynamic adjacency matrix implementation with:
  - File I/O support (`.txt`/`.dot`)
  - Validation for connected weighted graphs
  - Export to Graphviz format

- **Algorithm Library**  
  Implementations of:
  - DFS/BFS traversal (with custom Stack/Queue)
  - Shortest Path:
    - Dijkstra's Algorithm (single-pair)
    - Floyd-Warshall (all-pairs)
  - Minimum Spanning Tree (Prim's Algorithm)
  - Traveling Salesman Problem (Ant Colony Optimization)

### Advanced Features

- **Custom Data Structures**
  - Non-STL Stack/Queue implementations
  - Memory-efficient matrix operations

- **Console Interface**

  ```bash
  ==== Graph Manager ====
  Load graph from file...
  Enter filename: ./dot_examples/adj_matrix.txt 
  Graph loaded successfully.

  ==== Graph Manager ====
  1. Load new graph from file
  2. Traverse the graph in breadth and print the result
  3. Traverse the graph in depth and print the result
  4. Find the shortest path between any two vertices and print the result
  5. Find the shortest paths between all pairs of vertices in the graph and print the result matrix
  6. Search for the minimum spanning tree in the graph and print the resulting adjacency matrix
  7. Solve the Salesman problem, with output of the resulting route and its length
  8. Bonus. Comparison of methods for solving the traveling salesman problem
  0. Exit
  =======================
  Enter your choice:
  ```
