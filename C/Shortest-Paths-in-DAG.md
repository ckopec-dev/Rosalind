# Rosalind Problem: Shortest Paths in DAG

## Problem Understanding

We need to find the shortest paths from a source node to all other nodes in a directed acyclic graph (DAG). This is a classic single-source shortest path problem that can be solved efficiently in DAGs using topological sorting followed by dynamic programming.

## Solution Approach

1. **Topological Sort**: Since the graph is a DAG, we can use topological sorting to process vertices in the correct order
2. **Dynamic Programming**: Process vertices in topological order and update distances using relaxation
3. **Distance Initialization**: Initialize source distance as 0 and all others as infinity

## Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define MAX_NODES 1000
#define INF INT_MAX

// Structure to represent an edge
typedef struct Edge {
    int dest;
    int weight;
    struct Edge* next;
} Edge;

// Structure to represent a graph
typedef struct Graph {
    int num_vertices;
    Edge* adj_list[MAX_NODES];
} Graph;

// Function to create a new edge
Edge* create_edge(int dest, int weight) {
    Edge* edge = (Edge*)malloc(sizeof(Edge));
    edge->dest = dest;
    edge->weight = weight;
    edge->next = NULL;
    return edge;
}

// Function to add an edge to the graph
void add_edge(Graph* graph, int src, int dest, int weight) {
    Edge* edge = create_edge(dest, weight);
    edge->next = graph->adj_list[src];
    graph->adj_list[src] = edge;
}

// Function to perform topological sort using DFS
void topological_sort_util(Graph* graph, int vertex, int visited[], int stack[], int* stack_index) {
    visited[vertex] = 1;
    
    Edge* current = graph->adj_list[vertex];
    while (current != NULL) {
        if (!visited[current->dest]) {
            topological_sort_util(graph, current->dest, visited, stack, stack_index);
        }
        current = current->next;
    }
    
    stack[++(*stack_index)] = vertex;
}

// Function to get topological order
void topological_sort(Graph* graph, int order[]) {
    int visited[MAX_NODES] = {0};
    int stack[MAX_NODES];
    int stack_index = -1;
    
    for (int i = 0; i < graph->num_vertices; i++) {
        if (!visited[i]) {
            topological_sort_util(graph, i, visited, stack, &stack_index);
        }
    }
    
    // Copy the stack to order array in reverse order
    for (int i = 0; i <= stack_index; i++) {
        order[i] = stack[stack_index - i];
    }
}

// Function to find shortest paths from source in DAG
void shortest_paths_dag(Graph* graph, int source, int distances[]) {
    // Initialize distances
    for (int i = 0; i < graph->num_vertices; i++) {
        distances[i] = INF;
    }
    distances[source] = 0;
    
    // Get topological order
    int topo_order[MAX_NODES];
    topological_sort(graph, topo_order);
    
    // Process vertices in topological order
    for (int i = 0; i < graph->num_vertices; i++) {
        int u = topo_order[i];
        
        // If we've reached a vertex that's unreachable from source
        if (distances[u] == INF) continue;
        
        // Update distances of all adjacent vertices
        Edge* current = graph->adj_list[u];
        while (current != NULL) {
            int v = current->dest;
            int weight = current->weight;
            
            // Relax the edge
            if (distances[u] + weight < distances[v]) {
                distances[v] = distances[u] + weight;
            }
            
            current = current->next;
        }
    }
}

// Function to print distances
void print_distances(int distances[], int num_vertices) {
    for (int i = 0; i < num_vertices; i++) {
        if (distances[i] == INF) {
            printf("x ");
        } else {
            printf("%d ", distances[i]);
        }
    }
    printf("\n");
}

int main() {
    // Read input
    int num_vertices, num_edges;
    scanf("%d %d", &num_vertices, &num_edges);
    
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->num_vertices = num_vertices;
    
    // Initialize adjacency lists
    for (int i = 0; i < MAX_NODES; i++) {
        graph->adj_list[i] = NULL;
    }
    
    // Read edges
    for (int i = 0; i < num_edges; i++) {
        int src, dest, weight;
        scanf("%d %d %d", &src, &dest, &weight);
        add_edge(graph, src - 1, dest - 1, weight); // Convert to 0-indexed
    }
    
    // Find shortest paths from vertex 1 (0-indexed)
    int distances[MAX_NODES];
    shortest_paths_dag(graph, 0, distances);
    
    // Print result
    print_distances(distances, num_vertices);
    
    // Free memory
    for (int i = 0; i < graph->num_vertices; i++) {
        Edge* current = graph->adj_list[i];
        while (current != NULL) {
            Edge* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(graph);
    
    return 0;
}
```

## Key Points

1. **Topological Sorting**: Since we have a DAG, we can use DFS-based topological sort to process vertices in the correct order
2. **Efficiency**: This approach runs in O(V + E) time, which is optimal for DAGs
3. **Distance Initialization**: Source vertex gets distance 0, all others get infinity
4. **Relaxation**: For each vertex in topological order, we update distances to its neighbors

## Input Format

The input consists of:
- First line: number of vertices (V) and number of edges (E)
- Next E lines: each contains three integers representing an edge (source, destination, weight)

## Output Format

Space-separated shortest distances from source vertex to all other vertices. Use "x" for unreachable vertices.

## Example

For input:
```
5 6
1 2 10
1 3 3
2 4 2
3 2 4
3 4 8
3 5 2
```

Output would be the shortest distances from vertex 1 to all other vertices.