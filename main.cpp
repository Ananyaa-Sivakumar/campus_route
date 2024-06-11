#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <ctime>
#include <functional>

using namespace std;

const double INF = numeric_limits<double>::infinity();

struct Edge {
    string src;
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
        nodes[srcIndex].edges.push_back({src, dest, weight});
        nodes[destIndex].edges.push_back({dest, src, weight});
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
            mstEdges.insert({edge.src, edge.dest});
            mstEdges.insert({edge.dest, edge.src});
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

        set<pair<string, string>> addedEdges;
        for (const auto& node : nodes) {
            for (const auto& edge : node.edges) {
                if (addedEdges.find({node.name, edge.dest}) == addedEdges.end() &&
                    addedEdges.find({edge.dest, node.name}) == addedEdges.end()) {
                    file << "    \"" << node.name << "\" -- \"" << edge.dest << "\" [label=\"" << edge.weight << "\"";

                    if (pathEdges.find({node.name, edge.dest}) != pathEdges.end() ||
                        pathEdges.find({edge.dest, node.name}) != pathEdges.end()) {
                        file << ", color=blue, penwidth=2.0";
                    } else if (mstEdges.find({node.name, edge.dest}) != mstEdges.end() ||
                               mstEdges.find({edge.dest, node.name}) != mstEdges.end()) {
                        file << ", color=green, penwidth=2.0";
                    }

                    file << "];" << endl;
                    addedEdges.insert({node.name, edge.dest});
                }
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

    cout << "Updating traffic weights by adding random values to each edge...\n";
    
    // Iterate through all nodes and their edges
    for (auto& node : nodes) {
        for (auto& edge : node.edges) {
            double oldWeight = edge.weight;
            double randomValue = 1 + rand() % 10; // Random value between 1 and 10
            edge.weight += randomValue; // Add random value to the original weight
            cout << "Traffic weight from " << edge.src << " to " << edge.dest << " changed from " << oldWeight << " to " << edge.weight << " (added " << randomValue << ")\n";
            
            // Update the corresponding reverse edge in the destination node
            for (auto& revEdge : nodes[nodeIndices[edge.dest]].edges) {
                if (revEdge.dest == node.name) {
                    revEdge.weight = edge.weight;
                    break;
                }
            }
        }
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

    for (size_t i = 0; i < graph.getNodes().size(); ++i) {
        for (const auto& edge : graph.getNodes()[i].edges) {
            int u = nodeIndices.at(graph.getNodes()[i].name);
            int v = nodeIndices.at(edge.dest);
            if (distances[u] + edge.weight == distances[v]) {
                previous[v] = u;
            }
        }
    }

    int at = nodeIndices.at(end);
    while (at != -1) {
        path.push_back(graph.getNodes()[at].name);
        at = previous[at];
    }
    reverse(path.begin(), path.end());

    return path;
}

vector<string> shortestPathViaAnyNode(Graph& graph, const string& start, const vector<string>& viaNodes, const string& end, double& totalWeight) {
    vector<string> path;
    totalWeight = 0;

    path.push_back(start);

    for (size_t i = 0; i < viaNodes.size(); ++i) {
        double segmentWeight = 0;
        vector<string> segmentPath;

        if (i == 0) {
            segmentPath = shortestPath(graph, start, viaNodes[i], segmentWeight);
        } else {
            segmentPath = shortestPath(graph, viaNodes[i - 1], viaNodes[i], segmentWeight);
        }

        if (segmentPath.empty()) {
            cerr << "No path found from " << (i == 0 ? start : viaNodes[i - 1]) << " to " << viaNodes[i] << "." << endl;
            return {};
        }

        path.insert(path.end(), segmentPath.begin() + 1, segmentPath.end());
        totalWeight += segmentWeight;
    }

    double lastSegmentWeight = 0;
    vector<string> lastSegmentPath = shortestPath(graph, viaNodes.back(), end, lastSegmentWeight);
    if (lastSegmentPath.empty()) {
        cerr << "No path found from " << viaNodes.back() << " to " << end << "." << endl;
        return {};
    }

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

void printGraphWithWeights(const vector<Edge>& mst) {
    // Create an adjacency list to store the MST
    unordered_map<string, vector<pair<string, double>>> adjList;

    for (const auto& edge : mst) {
        adjList[edge.src].push_back({edge.dest, edge.weight});
        adjList[edge.dest].push_back({edge.src, edge.weight});
    }

    // Start from any node (we choose the first node in the mst list)
    string startNode = mst.front().src;

    // Use a DFS to print the MST edges
    set<string> visited;
    function<void(const string&)> dfs = [&](const string& node) {
        visited.insert(node);
        cout << node << " -> ";
        for (const auto& neighbor : adjList[node]) {
            if (visited.find(neighbor.first) == visited.end()) {
                cout << neighbor.first << " (" << neighbor.second << "), ";
                dfs(neighbor.first);
            }
        }
    };

    dfs(startNode);
    cout << "END" << endl;
}

vector<Edge> calculateMST(Graph& graph, const string& start, double& totalWeight) {
    vector<Edge> mst = graph.primsMST(start);
    totalWeight = 0;
    for (auto& edge : mst) {
        totalWeight += edge.weight;
    }
    return mst;
}

int main() {
    Graph campusGraph;
    loadGraphFromFile(campusGraph, "graph_data.txt");

    int choice;
    do {
        cout << "\n========== Campus Navigation System ==========\n";
        cout << "1) Shortest Path\n";
        cout << "2) Optimal Route for Campus Visit (MST)\n";
        cout << "3) Display Graph\n";
        cout << "4) Visualize Graph\n";
        cout << "5) Update Traffic Weights\n";
        cout << "0) Exit\n";
        cout << "==============================================\n";
        cout << "Enter your choice: ";
        cin >> choice;
        cout << endl;

        if (choice == 1) {
            int subChoice;
            cout << "========== Shortest Path Options ==========\n";
            cout << "1) Normal Shortest Path\n";
            cout << "2) Shortest Path via Any Node\n";
            cout << "===========================================\n";
            cout << "Enter your choice: ";
            cin >> subChoice;
            cout << endl;

            string start, end, via;
            double totalWeight;
            vector<string> path;

            if (subChoice == 1) {
                cout << "Enter start node: ";
                cin >> start;
                cout << "Enter end node: ";
                cin >> end;
                cout << endl;

                if (!campusGraph.hasNode(start) || !campusGraph.hasNode(end)) {
                    cerr << "Error: Start or end node not found in the graph." << endl;
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
                cout << endl;

                if (!campusGraph.hasNode(start) || !campusGraph.hasNode(end)) {
                    cerr << "Error: Start or end node not found in the graph." << endl;
                    continue;
                }

                path = shortestPathViaAnyNode(campusGraph, start, {via}, end, totalWeight);
            }

            if (!path.empty()) {
                cout << "========== Shortest Path Result ==========\n";
                cout << "Path: ";
                for (const auto& node : path) {
                    cout << node << " ";
                }
                cout << "\nTotal weight: " << totalWeight << endl;
                cout << "===========================================\n";
                campusGraph.visualizeGraph("graph.dot", path, {});
            }

        } else if (choice == 2) {
            string start;
            cout << "Enter starting node: ";
            cin >> start;
            cout << endl;

            try {
                double totalWeight = 0;
                vector<Edge> mst = calculateMST(campusGraph, start, totalWeight);
                cout << "========== Minimum Spanning Tree ==========\n";
                printGraphWithWeights(mst);
                cout << "Total weight of MST: " << totalWeight << endl;
                cout << "===========================================\n";
                campusGraph.visualizeGraph("graph.dot", {}, mst);

            } catch (const out_of_range& e) {
                cerr << "Error: " << e.what() << endl;
            }

        } else if (choice == 3) {
            cout << "========== Graph Structure ==========\n";
            campusGraph.displayGraph();
            cout << "=====================================\n";

        } else if (choice == 4) {
            cout << "Generating graph visualization...\n";
            campusGraph.visualizeGraph("graph.dot");

        } else if (choice == 5) {
            campusGraph.updateTrafficWeights();
            cout << "Traffic weights updated randomly.\n";

         
        } else if (choice == 0) {
            cout << "Exiting...\n";
        } else {
            cout << "Invalid choice. Please enter a number between 0 and 5.\n";
        }

    } while (choice != 0);

    return 0;
}



