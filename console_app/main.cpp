#include <iostream>
#include <s21_console_app.hpp>

int main() {
  s21::ConsoleApp conApp;
  cout << "\n==== Graph Manager ====\n";
  cout << "Load graph from file...\n";
  while (!conApp.SetFileNameFromConsole()) {
    cout << "\n==== Graph Manager ====\n";
    cout << "Load graph from file...\n";
  }
  int vertex1, vertex2, startVertex, N;

  while (true) {
    conApp.DisplayMenu();
    int choice;
    cin >> choice;

    switch (choice) {
      case 1:
        conApp.SetFileNameFromConsole();
        break;
      case 2:
        cout << "Enter start vertex: ";
        cin >> startVertex;
        conApp.SetStartVertex(startVertex);
        conApp.ConsoleBFS();
        cin.ignore();
        cin.get();
        break;
      case 3:
        cout << "Enter start vertex: ";
        cin >> startVertex;
        conApp.SetStartVertex(startVertex);
        conApp.ConsoleDFS();
        cin.ignore();
        cin.get();
        break;
      case 4:
        cout << "Enter Vertex1 and Vertex2: ";
        cin >> vertex1 >> vertex2;
        conApp.SetVertices(vertex1, vertex2);
        conApp.ConsoleDijkstra();
        cin.ignore();
        cin.get();
        break;
      case 5:
        conApp.ConsoleFloydWarshall();
        cin.ignore();
        cin.get();
        break;
      case 6:
        conApp.ConsolePrim();
        cin.ignore();
        cin.get();
        break;
      case 7:
        conApp.ConsoleAntAlgo();
        cin.ignore();
        cin.get();
        break;
      case 8:
        cout << "Enter N (number of times in a row to solve the salesman's "
                "problem): ";
        cin >> N;
        conApp.SetN(N);
        conApp.ComparisonTests();
        cin.ignore();
        cin.get();
        break;
      case 0:
        std::cout << "Exiting...\n";
        return 0;
      default:
        std::cout << "Invalid choice. Please try again.\n";
        break;
    }
  }

  return 0;
}
