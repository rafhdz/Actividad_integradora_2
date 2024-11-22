#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <cmath>
#include <limits>

using namespace std;

// Estructura para las aristas
struct Edge {
    int u, v;
    double weight;
    Edge(int u, int v, double w) : u(u), v(v), weight(w) {}
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

// Estructura para los puntos (centrales)
struct Point {
    double x, y;
    Point() : x(0), y(0) {}  // Constructor por defecto
    Point(double x, double y) : x(x), y(y) {}
};

// Disjoint Set Union (Union-Find) para Kruskal
class DSU {
    vector<int> parent, rank;
public:
    DSU(int n) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int find(int x) {
        if (parent[x] != x)
            parent[x] = find(parent[x]);
        return parent[x];
    }
    bool unite(int x, int y) {
        int xr = find(x), yr = find(y);
        if (xr == yr) return false;
        if (rank[xr] < rank[yr]) parent[xr] = yr;
        else if (rank[xr] > rank[yr]) parent[yr] = xr;
        else { parent[yr] = xr; rank[xr]++; }
        return true;
    }
};

// Estructura para almacenar el grafo de flujo
struct FlowEdge {
    int v, rev;
    int capacity, flow;
    FlowEdge(int v, int rev, int capacity) : v(v), rev(rev), capacity(capacity), flow(0) {}
};

class MaxFlow {
    int N;
    vector<vector<FlowEdge> > adj;
    vector<int> level, ptr;
public:
    MaxFlow(int N) : N(N), adj(N), level(N), ptr(N) {}

    void addEdge(int u, int v, int capacity) {
        adj[u].emplace_back(v, adj[v].size(), capacity);
        adj[v].emplace_back(u, adj[u].size() - 1, 0);
    }

    bool bfs(int s, int t) {
        fill(level.begin(), level.end(), -1);
        level[s] = 0;
        queue<int> q;
        q.push(s);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (const FlowEdge& e : adj[u]) {
                if (level[e.v] == -1 && e.flow < e.capacity) {
                    level[e.v] = level[u] + 1;
                    q.push(e.v);
                }
            }
        }
        return level[t] != -1;
    }

    int dfs(int u, int t, int pushed) {
        if (pushed == 0) return 0;
        if (u == t) return pushed;

        for (int& cid = ptr[u]; cid < adj[u].size(); ++cid) {
            FlowEdge& e = adj[u][cid];
            if (level[u] + 1 != level[e.v] || e.flow == e.capacity) continue;
            int tr = dfs(e.v, t, min(pushed, e.capacity - e.flow));
            if (tr == 0) continue;
            e.flow += tr;
            adj[e.v][e.rev].flow -= tr;
            return tr;
        }
        return 0;
    }

    int maxFlow(int s, int t) {
        int flow = 0;
        while (bfs(s, t)) {
            fill(ptr.begin(), ptr.end(), 0);
            while (int pushed = dfs(s, t, numeric_limits<int>::max())) {
                flow += pushed;
            }
        }
        return flow;
    }
};

// Algoritmo 2-opt para mejorar la ruta del TSP
void twoOpt(vector<int>& route, const vector<vector<double> >& distance) {
    int N = route.size() - 1; // route[0] == route[N]
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 1; i < N - 1; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double delta = (distance[route[i - 1]][route[j]] + distance[route[i]][route[(j + 1) % N]]) -
                               (distance[route[i - 1]][route[i]] + distance[route[j]][route[(j + 1) % N]]);
                if (delta < -1e-6) {
                    reverse(route.begin() + i, route.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }
}

// Algoritmo de Vecino Más Cercano para el TSP
vector<int> nearestNeighborTSP(const vector<vector<double> >& distance) {
    int N = distance.size();
    vector<bool> visited(N, false);
    vector<int> route;
    route.push_back(0); // Comenzar desde el nodo 0
    visited[0] = true;

    for (int i = 1; i < N; ++i) {
        int current = route.back();
        int next = -1;
        double minDist = numeric_limits<double>::max();

        for (int j = 0; j < N; ++j) {
            if (!visited[j] && distance[current][j] < minDist) {
                minDist = distance[current][j];
                next = j;
            }
        }
        if (next != -1) {
            route.push_back(next);
            visited[next] = true;
        }
    }
    route.push_back(0); // Regresar al inicio
    return route;
}

// Función para calcular las celdas de Voronoi (simplificado)
vector<vector<Point> > computeVoronoiCells(const vector<Point>& points) {
    vector<vector<Point> > cells(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        cells[i].push_back(Point(points[i].x - 50, points[i].y - 50));
        cells[i].push_back(Point(points[i].x - 50, points[i].y + 50));
        cells[i].push_back(Point(points[i].x + 50, points[i].y + 50));
        cells[i].push_back(Point(points[i].x + 50, points[i].y - 50));
    }
    return cells;
}

int main() {
    ifstream inputFile("input.txt");
    if (!inputFile.is_open()) {
        cerr << "No se pudo abrir el archivo de entrada.\n";
        return 1;
    }

    int N;
    inputFile >> N;

    vector<vector<double> > distance(N, vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            inputFile >> distance[i][j];

    vector<vector<int> > capacity(N, vector<int>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            inputFile >> capacity[i][j];

    vector<Point> centrals(N);
    for (int i = 0; i < N; ++i)
        inputFile >> centrals[i].x >> centrals[i].y;

    inputFile.close();

    // Parte 1: Árbol de Expansión Mínima
    vector<Edge> edges;
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            edges.emplace_back(i, j, distance[i][j]);

    sort(edges.begin(), edges.end());
    DSU dsu(N);
    vector<pair<int, int> > mstEdges;
    for (const Edge& e : edges) {
        if (dsu.unite(e.u, e.v)) {
            mstEdges.emplace_back(e.u, e.v);
        }
    }

    cout << "1. Forma de cablear las colonias con fibra:\n";
    for (const auto& e : mstEdges) {
        char colony1 = 'A' + e.first;
        char colony2 = 'A' + e.second;
        cout << "(" << colony1 << "," << colony2 << ")\n";
    }

    // Parte 2: Ruta del Viajante con mejora 2-opt
    vector<int> route = nearestNeighborTSP(distance);
    twoOpt(route, distance);
    cout << "2. Ruta de correspondencia:\n";
    for (size_t i = 0; i < route.size(); ++i) {
        if (i > 0) cout << "-";
        cout << char('A' + route[i]);
    }
    cout << "\n";

    // Parte 3: Flujo Máximo
    MaxFlow mf(N);
    for (int u = 0; u < N; ++u)
        for (int v = 0; v < N; ++v)
            if (capacity[u][v] > 0)
                mf.addEdge(u, v, capacity[u][v]);

    cout << "3. Flujo máximo: " << mf.maxFlow(0, N - 1) << "\n";

    // Parte 4: Voronoi (simplificado)
    auto voronoiCells = computeVoronoiCells(centrals);
    cout << "4. Polígonos de Voronoi:\n";
    for (size_t i = 0; i < voronoiCells.size(); ++i) {
        cout << "Polígono " << i + 1 << ":\n";
        for (const Point& p : voronoiCells[i]) {
            cout << "(" << p.x << "," << p.y << ") ";
        }
        cout << "\n";
    }

    return 0;
}
