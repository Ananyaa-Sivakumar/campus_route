#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <cstdlib>  // For system(), rand(), srand()
#include <ctime>    // For time()

using namespace std;

const double INF = numeric_limits<double>::infinity();

struct Edge {
    string dest;
    double weight;

    bool operator>(const Edge& other) const {
        return weight > other.weight;
    }
};

struct Node {
    string name;
    vector<Edge> edges;
};

class Graph {
private:
    vector<Node> nodes;
    unordered_map<string, int> nodeIndices;

public:
    void addEdge(const string& src, const string& dest, double weight) {
        int srcIndex = nodeIndices[src];
        int destIndex = nodeIndices[dest];
        nodes[srcIndex].edges.push_back({dest, weight});
        nodes[destIndex].edges.push_back({src, weight});
    }

    void addNode(const string& name) {
        if (nodeIndices.find(name) == nodeIndices.end()) {
            nodeIndices[name] = nodes.size();
            nodes.push_back({name, {}});
        }
    }

    bool hasNode(const string& name) const {
        return nodeIndices.find(name) != nodeIndices.end();
    }

    const vector<Node>& getNodes() const {
        return nodes;
    }

    const unordered_map<string, int>& getNodeIndices() const {
        return nodeIndices;
    }

    vector<double> dijkstra(const string& start) {
        vector<double> distances(nodes.size(), INF);
        distances[nodeIndices[start]] = 0;

        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
        pq.push(make_pair(0, nodeIndices[start]));

        while (!pq.empty()) {
            double dist = pq.top().first;
            int current = pq.top().second;
            pq.pop();

            if (dist > distances[current]) continue;

            for (const auto& edge : nodes[current].edges) {
                double newDist = dist + edge.weight;
                if (newDist < distances[nodeIndices[edge.dest]]) {
                    distances[nodeIndices[edge.dest]] = newDist;
                    pq.push(make_pair(newDist, nodeIndices[edge.dest]));
                }
            }
        }

        return distances;
    }

    vector<Edge> primsMST(const string& start) {
        vector<bool> inMST(nodes.size(), false);
        vector<Edge> result;

        priority_queue<Edge, vector<Edge>, greater<Edge>> pq;

        for (const auto& edge : nodes[nodeIndices[start]].edges) {
            pq.push(edge);
        }
        inMST[nodeIndices[start]] = true;

        while (!pq.empty()) {
            Edge edge = pq.top();
            pq.pop();

            if (inMST[nodeIndices[edge.dest]]) continue;

            inMST[nodeIndices[edge.dest]] = true;
            result.push_back(edge);

            for (const auto& nextEdge : nodes[nodeIndices[edge.dest]].edges) {
                if (!inMST[nodeIndices[nextEdge.dest]]) {
                    pq.push(nextEdge);
                }
            }
        }

        return result;
    }

    void generateDOTFile(const string& filename, const vector<string>& path = {}, const vector<Edge>& mst = {}, const string& startNode = "") const {
        ofstream file(filename);
        file << "graph G {" << endl;

        set<pair<string, string>> mstEdges;
        for (const auto& edge : mst) {
            mstEdges.insert({edge.dest, edge.dest});
        }

        set<pair<string, string>> pathEdges;
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            pathEdges.insert({path[i], path[i + 1]});
            pathEdges.insert({path[i + 1], path[i]});
        }

        for (const auto& node : nodes) {
            file << "    \"" << node.name << "\"";
            if (node.name == startNode) {
                file << " [color=red]";
            }
            file << ";" << endl;
        }

        for (const auto& node : nodes) {
            for (const auto& edge : node.edges) {
                file << "    \"" << node.name << "\" -- \"" << edge.dest << "\" [label=\"" << edge.weight << "\"";

                if (pathEdges.find({node.name, edge.dest}) != pathEdges.end()) {
                    file << ", color=blue, penwidth=2.0";
                } else if (mstEdges.find({node.name, edge.dest}) != mstEdges.end()) {
                    file << ", color=green, penwidth=2.0";
                }

                file << "];" << endl;
            }
        }
        file << "}" << endl;
        file.close();
    }

    void visualizeGraph(const string& filename, const vector<string>& path = {}, const vector<Edge>& mst = {}) const {
        generateDOTFile(filename, path, mst);
        string command = "dot -Tpng " + filename + " -o graph.png && start graph.png";
        system(command.c_str());
    }

    void displayGraph() const {
        for (const auto& node : nodes) {
            cout << node.name << " -> ";
            for (const auto& edge : node.edges) {
                cout << "(" << edge.dest << ", " << edge.weight << ") ";
            }
            cout << endl;
        }
    }

    
void updateTrafficWeights() {
    srand(time(0)); // Initialize random seed

    // Select a random node
    int randomNodeIndex = rand() % nodes.size();
    Node& randomNode = nodes[randomNodeIndex];

    cout << "Randomly selected node: " << randomNode.name << endl;

    // Update traffic weights for the edges of the selected node
    for (auto& edge : randomNode.edges) {
        double oldWeight = edge.weight;
        edge.weight = 1 + rand() % 10; // Random weight between 1 and 10
        cout << "Traffic weight from " << randomNode.name << " to " << edge.dest << " changed from " << oldWeight << " to " << edge.weight << endl;
    }
}


};

vector<string> shortestPath(Graph& graph, const string& start, const string& end, double& totalWeight) {
    vector<double> distances = graph.dijkstra(start);
    vector<int> previous(graph.getNodes().size(), -1);
    const unordered_map<string, int>& nodeIndices = graph.getNodeIndices();

    if (distances[nodeIndices.at(end)] == INF) {
        cerr << "No path found from " << start << " to " << end << "." << endl;
        return {};
    }

    vector<string> path;
    totalWeight = distances[nodeIndices.at(end)];

    // Update the 'previous' array during Dijkstra's algorithm execution
    for (size_t i = 0; i < graph.getNodes().size(); ++i) {
        for (const auto& edge : graph.getNodes()[i].edges) {
            int u = nodeIndices.at(graph.getNodes()[i].name);
            int v = nodeIndices.at(edge.dest);
            if (distances[u] + edge.weight == distances[v]) {
                previous[v] = u;
            }
        }
    }

    // Reconstruct the shortest path using the 'previous' array
    int at = nodeIndices.at(end);
    while (at != -1) {
        path.push_back(graph.getNodes()[at].name);
        at = previous[at];
    }
    reverse(path.begin(), path.end());

    return path;
}


vector<string> shortestPathViaAnyNode(Graph& graph, const string& start, const vector<string>& viaNodes, const string& end, double& totalWeight) {
    // Initialize path and total weight
    vector<string> path;
    totalWeight = 0;

    // Start with the start node
    path.push_back(start);

    // Iterate through via nodes and find shortest path between them
    for (size_t i = 0; i < viaNodes.size(); ++i) {
        double segmentWeight = 0;
        vector<string> segmentPath;

        if (i == 0) {
            // Find shortest path from start to first via node
            segmentPath = shortestPath(graph, start, viaNodes[i], segmentWeight);
        } else {
            // Find shortest path between via nodes
            segmentPath = shortestPath(graph, viaNodes[i - 1], viaNodes[i], segmentWeight);
        }

        if (segmentPath.empty()) {
            cerr << "No path found from " << (i == 0 ? start : viaNodes[i - 1]) << " to " << viaNodes[i] << "." << endl;
            return {};
        }

        // Append segment path to total path
        path.insert(path.end(), segmentPath.begin() + 1, segmentPath.end());
        totalWeight += segmentWeight;
    }

    // Find shortest path from last via node to end node
    double lastSegmentWeight = 0;
    vector<string> lastSegmentPath = shortestPath(graph, viaNodes.back(), end, lastSegmentWeight);
    if (lastSegmentPath.empty()) {
        cerr << "No path found from " << viaNodes.back() << " to " << end << "." << endl;
        return {};
    }

    // Append last segment path to total path
    path.insert(path.end(), lastSegmentPath.begin() + 1, lastSegmentPath.end());
    totalWeight += lastSegmentWeight;

    return path;
}

void loadGraphFromFile(Graph& graph, const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    int numEdges;
    file >> numEdges;

    for (int i = 0; i < numEdges; ++i) {
        string src, dest;
        double weight;
        file >> src >> dest >> weight;
        graph.addNode(src);
        graph.addNode(dest);
        graph.addEdge(src, dest, weight);
    }

    file.close();
}

void displayNodes(const Graph& graph) {
    cout << "Nodes in the graph: ";
    for (const auto& node : graph.getNodes()) {
        cout << node.name << " ";
    }
    cout << endl;
}

void printGraphWithWeights(const vector<Edge>& mst, const string& start) {
    cout << "Optimal Route for Entire Campus (MST):\n";
    for (const auto& edge : mst) {
        cout << start << " -> " << edge.dest << " [weight: " << edge.weight << "]\n";
    }
}

vector<Edge> calculateMST(Graph& graph, const string& start, double& totalWeight) {
    vector<Edge> mst = graph.primsMST(start);
    totalWeight = 0;
    for (const auto& edge : mst) {
        totalWeight += edge.weight;
    }
    return mst;
}

int main() {
    Graph campusGraph;
    loadGraphFromFile(campusGraph, "graph_data.txt");

    int choice;
    do {
        cout << "Menu:\n";
        cout << "1) Shortest Path\n";
        cout << "2) Optimal Route for Campus Visit (MST)\n";
        cout << "3) Display Graph\n";
        cout << "4) Visualize Graph\n";
        cout << "5) Update Traffic Weights\n"; // New menu option
        cout << "0) Exit\n";
        cout << "Enter your choice: ";
        cin >> choice;

        if (choice == 1) {
            int subChoice;
            cout << "1) Normal Shortest Path\n";
            cout << "2) Shortest Path via Any Node\n";
            cout << "Enter your choice: ";
            cin >> subChoice;

            string start, end, via;
            double totalWeight;
            vector<string> path;

            if (subChoice == 1) {
                cout << "Enter start node: ";
                cin >> start;
                cout << "Enter end node: ";
                cin >> end;
                if (!campusGraph.hasNode(start) || !campusGraph.hasNode(end)) {
                    cerr << "Start or end node not found in the graph." << endl;
                    continue;
                }

                path = shortestPath(campusGraph, start, end, totalWeight);
            } else if (subChoice == 2) {
                cout << "Enter start node: ";
                cin >> start;
                cout << "Enter via node: ";
                cin >> via;
                cout << "Enter end node: ";
                cin >> end;
                if (!campusGraph.hasNode(start) || !campusGraph.hasNode(end)) {
                    cerr << "Start or end node not found in the graph." << endl;
                    continue;
                }

                path = shortestPathViaAnyNode(campusGraph, start, {via}, end, totalWeight);
            }

            if (!path.empty()) {
                cout << "Shortest path: ";
                for (const auto& node : path) {
                    cout << node << " ";
                }
                cout << "\nTotal weight: " << totalWeight << endl;
                campusGraph.visualizeGraph("graph.dot", path, {});
            }

        } else if (choice == 2) {
            // Optimal Route for Campus Visit (MST) option
            string start;
            cout << "Enter starting node: ";
            cin >> start;

            try {
                double totalWeight = 0;
                vector<Edge> mst = calculateMST(campusGraph, start, totalWeight);
                printGraphWithWeights(mst, start);
                cout << "Total weight of MST: " << totalWeight << endl;
                campusGraph.visualizeGraph("graph.dot", {}, mst);

            } catch (const out_of_range& e) {
                cerr << e.what() << endl;
            }
        } else if (choice == 3) {
            campusGraph.displayGraph();

        } else if (choice == 4) {
            campusGraph.visualizeGraph("graph.dot");

        } else if (choice == 5) {
            // Update Traffic Weights option
            campusGraph.updateTrafficWeights();
            cout << "Traffic weights updated randomly.\n";

        } else if (choice == 0) {
            cout << "Exiting..." << endl;
        } else {
            cout<< "Invalid choice. Please try again." << endl;
        }

    } while (choice != 0);

    return 0;
}

