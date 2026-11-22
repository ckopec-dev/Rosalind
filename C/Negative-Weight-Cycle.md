# Rosalind Problem: Negative Weight Cycle (C Solution)

## Problem Understanding

The negative weight cycle problem asks us to detect if there's a negative weight cycle in a weighted directed graph. A negative weight cycle is a cycle where the sum of edge weights is negative, which can cause issues in shortest path algorithms.

## Approach

I'll use the Bellman-Ford algorithm with modifications to detect negative cycles:
1. Run Bellman-Ford for V-1 iterations
2. Run one more iteration to check if any distances can still be relaxed
3. If they can, a negative cycle exists

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define MAX_NODES 1000
#define INF INT_MAX

// Edge structure
typedef struct {
    int source;
    int destination;
    int weight;
} Edge;

// Function to detect negative weight cycle using modified Bellman-Ford
int hasNegativeCycle(Edge edges[], int numEdges, int numNodes) {
    // Initialize distances
    int dist[MAX_NODES];
    for (int i = 0; i < numNodes; i++) {
        dist[i] = INF;
    }
    dist[0] = 0; // Start from node 0
    
    // Relax edges V-1 times
    for (int i = 0; i < numNodes - 1; i++) {
        for (int j = 0; j < numEdges; j++) {
            int u = edges[j].source;
            int v = edges[j].destination;
            int weight = edges[j].weight;
            
            if (dist[u] != INF && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
            }
        }
    }
    
    // Check for negative cycles - run one more iteration
    for (int j = 0; j < numEdges; j++) {
        int u = edges[j].source;
        int v = edges[j].destination;
        int weight = edges[j].weight;
        
        if (dist[u] != INF && dist[u] + weight < dist[v]) {
            return 1; // Negative cycle detected
        }
    }
    
    return 0; // No negative cycle
}

int main() {
    int numNodes, numEdges;
    
    // Read input
    scanf("%d %d", &numNodes, &numEdges);
    
    Edge *edges = (Edge*)malloc(numEdges * sizeof(Edge));
    
    // Read edges
    for (int i = 0; i < numEdges; i++) {
        int u, v, w;
        scanf("%d %d %d", &u, &v, &w);
        edges[i].source = u - 1;      // Convert to 0-based indexing
        edges[i].destination = v - 1;
        edges[i].weight = w;
    }
    
    // Check for negative cycle
    if (hasNegativeCycle(edges, numEdges, numNodes)) {
        printf("1\n");  // Negative cycle exists
    } else {
        printf("0\n");  // No negative cycle
    }
    
    free(edges);
    return 0;
}
```

## Explanation

### Algorithm Steps:
1. **Initialization**: Set all distances to infinity except the source (set to 0)
2. **Relaxation**: Run Bellman-Ford for V-1 iterations to find shortest paths
3. **Negative Cycle Detection**: Run one more iteration to check if any distances can still be reduced
4. **Result**: If any distance is updated in the final iteration, a negative cycle exists

### Key Points:
- Uses 0-based indexing internally but reads 1-based input
- Time complexity: O(VE) where V is vertices and E is edges
- Space complexity: O(V + E)
- The algorithm works because if there's a negative cycle, we can always find a shorter path by going around the cycle

### Input Format:
```
Number_of_nodes Number_of_edges
Source Destination Weight
...
```

### Output:
- `1` if negative weight cycle exists
- `0` if no negative weight cycle exists

This solution handles the standard Bellman-Ford algorithm modification needed to detect negative cycles in directed graphs.