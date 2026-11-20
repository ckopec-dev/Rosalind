# Rosalind Problem: Square_in_a_Graph

## Problem Understanding

In graph theory, a square in a graph is a cycle of length 4 (a 4-cycle). Given an undirected graph, we need to determine if there exists a 4-cycle in the graph.

## Approach

I'll use an adjacency matrix representation and check for 4-cycles by:
1. For each pair of vertices that are connected by an edge
2. Check if they share at least two common neighbors
3. If so, we have found a 4-cycle

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_N 100

int adj_matrix[MAX_N][MAX_N];
int n, m;

// Function to check if there's a 4-cycle in the graph
int has_square() {
    // For each edge (u,v) in the graph
    for (int u = 0; u < n; u++) {
        for (int v = u + 1; v < n; v++) {
            // If u and v are connected
            if (adj_matrix[u][v]) {
                // Count common neighbors of u and v
                int common_neighbors = 0;
                int common_vertices[MAX_N];
                int idx = 0;
                
                // Find all vertices that are neighbors of both u and v
                for (int w = 0; w < n; w++) {
                    if (adj_matrix[u][w] && adj_matrix[v][w]) {
                        common_vertices[idx++] = w;
                        common_neighbors++;
                    }
                }
                
                // If there are at least 2 common neighbors, we have a 4-cycle
                // The 4-cycle would be: u - v - w1 - w2 - u (where w1 and w2 are 
                // two common neighbors)
                if (common_neighbors >= 2) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

int main() {
    // Read number of vertices and edges
    scanf("%d %d", &n, &m);
    
    // Initialize adjacency matrix
    memset(adj_matrix, 0, sizeof(adj_matrix));
    
    // Read edges
    for (int i = 0; i < m; i++) {
        int u, v;
        scanf("%d %d", &u, &v);
        // Convert to 0-indexed
        u--; v--;
        adj_matrix[u][v] = 1;
        adj_matrix[v][u] = 1;
    }
    
    // Check if there's a square (4-cycle)
    if (has_square()) {
        printf("1\n");
    } else {
        printf("0\n");
    }
    
    return 0;
}
```

## Explanation

The algorithm works as follows:

1. **Input Processing**: Read the number of vertices `n` and edges `m`, then read all edges and store them in an adjacency matrix.

2. **4-Cycle Detection**: For each edge (u,v) in the graph:
   - Find all vertices that are neighbors of both u and v (common neighbors)
   - If there are at least 2 common neighbors, we can form a 4-cycle
   - A 4-cycle would be: u → v → w1 → w2 → u where w1 and w2 are two common neighbors

3. **Output**: Return 1 if a 4-cycle exists, 0 otherwise.

## Time Complexity
- O(n³) where n is the number of vertices
- For each edge, we check all vertices to find common neighbors

## Space Complexity
- O(n²) for storing the adjacency matrix

This solution efficiently detects whether any 4-cycle exists in the given graph.