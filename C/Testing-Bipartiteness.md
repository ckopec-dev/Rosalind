# Rosalind Problem: Testing Bipartiteness (C Solution)

## Problem Understanding

A graph is bipartite if its vertices can be divided into two disjoint sets such that every edge connects a vertex from one set to a vertex from the other set. In other words, we can color the graph using only 2 colors such that no two adjacent vertices have the same color.

## Approach

I'll use BFS (Breadth-First Search) to check if the graph is bipartite:
1. Color each vertex with one of two colors (0 or 1)
2. Start from any unvisited vertex and color it with color 0
3. For each neighbor, if unvisited, color it with opposite color
4. If a neighbor is already visited and has the same color, the graph is not bipartite

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define MAX_VERTICES 1000

// Structure for adjacency list
typedef struct Node {
    int vertex;
    struct Node* next;
} Node;

// Structure for graph
typedef struct {
    int num_vertices;
    int num_edges;
    Node* adj_list[MAX_VERTICES];
} Graph;

// Function to create a new node
Node* create_node(int vertex) {
    Node* new_node = (Node*)malloc(sizeof(Node));
    new_node->vertex = vertex;
    new_node->next = NULL;
    return new_node;
}

// Function to create a graph
Graph* create_graph(int num_vertices) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->num_vertices = num_vertices;
    graph->num_edges = 0;
    
    for (int i = 0; i < MAX_VERTICES; i++) {
        graph->adj_list[i] = NULL;
    }
    
    return graph;
}

// Function to add edge to the graph
void add_edge(Graph* graph, int src, int dest) {
    // Add edge from src to dest
    Node* new_node = create_node(dest);
    new_node->next = graph->adj_list[src];
    graph->adj_list[src] = new_node;
    
    // Add edge from dest to src (for undirected graph)
    new_node = create_node(src);
    new_node->next = graph->adj_list[dest];
    graph->adj_list[dest] = new_node;
    
    graph->num_edges++;
}

// Function to check if graph is bipartite
bool is_bipartite(Graph* graph) {
    // Color array: -1 = unvisited, 0 = color 1, 1 = color 2
    int color[MAX_VERTICES];
    
    // Initialize all vertices as unvisited
    for (int i = 0; i < graph->num_vertices; i++) {
        color[i] = -1;
    }
    
    // BFS queue
    int queue[MAX_VERTICES];
    int front = 0, rear = 0;
    
    // Check all components of the graph
    for (int start_vertex = 0; start_vertex < graph->num_vertices; start_vertex++) {
        if (color[start_vertex] == -1) {
            // Start BFS from this vertex
            color[start_vertex] = 0;
            queue[rear++] = start_vertex;
            
            while (front < rear) {
                int current = queue[front++];
                
                // Visit all neighbors
                Node* neighbor = graph->adj_list[current];
                while (neighbor != NULL) {
                    int adj_vertex = neighbor->vertex;
                    
                    if (color[adj_vertex] == -1) {
                        // Color with opposite color
                        color[adj_vertex] = 1 - color[current];
                        queue[rear++] = adj_vertex;
                    } else if (color[adj_vertex] == color[current]) {
                        // Same color as current vertex - not bipartite
                        return false;
                    }
                    
                    neighbor = neighbor->next;
                }
            }
        }
    }
    
    return true;
}

int main() {
    int num_vertices, num_edges;
    
    // Read the number of vertices and edges
    scanf("%d %d", &num_vertices, &num_edges);
    
    // Create graph
    Graph* graph = create_graph(num_vertices);
    
    // Read edges
    for (int i = 0; i < num_edges; i++) {
        int src, dest;
        scanf("%d %d", &src, &dest);
        // Convert to 0-based indexing
        add_edge(graph, src - 1, dest - 1);
    }
    
    // Check if bipartite
    if (is_bipartite(graph)) {
        printf("1\n");  // Yes, bipartite
    } else {
        printf("-1\n"); // No, not bipartite
    }
    
    // Free memory
    for (int i = 0; i < graph->num_vertices; i++) {
        Node* current = graph->adj_list[i];
        while (current != NULL) {
            Node* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(graph);
    
    return 0;
}
```

## Explanation

1. **Data Structures**:
   - `Node`: Represents an adjacency list node
   - `Graph`: Contains the adjacency list representation of the graph
   - Color array: Tracks colors of vertices (-1 = unvisited, 0 or 1 = colors)

2. **Algorithm**:
   - Uses BFS traversal to color vertices
   - Start with any unvisited vertex and color it with color 0
   - Color all neighbors with opposite color (1 - current_color)
   - If we encounter a neighbor that's already colored with same color, graph is not bipartite

3. **Key Points**:
   - Handles disconnected components by checking all vertices
   - Uses 0-based indexing internally but reads 1-based input
   - Returns 1 for bipartite graphs, -1 for non-bipartite

## Time and Space Complexity

- **Time Complexity**: O(V + E) where V is vertices and E is edges
- **Space Complexity**: O(V + E) for adjacency list and color array

This solution efficiently determines whether a given undirected graph is bipartite using BFS traversal.