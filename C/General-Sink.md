# Rosalind Problem: General Sink (GEN_SINK)

## Problem Description

A sink is a node with no outgoing edges. In a directed graph, we want to find all nodes that are sinks.

## Solution Approach

1. Read the number of nodes and edges
2. Create an adjacency list representation of the graph
3. Identify nodes with no outgoing edges (sinks)
4. Output the sink nodes

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    
    // Create adjacency list - we'll use a boolean matrix to track edges
    int** graph = (int**)malloc((n + 1) * sizeof(int*));
    for (int i = 0; i <= n; i++) {
        graph[i] = (int*)calloc((n + 1), sizeof(int));
    }
    
    // Read edges and populate the graph
    for (int i = 0; i < m; i++) {
        int u, v;
        scanf("%d %d", &u, &v);
        graph[u][v] = 1;
    }
    
    // Find all sinks (nodes with no outgoing edges)
    int count = 0;
    int* sinks = (int*)malloc((n + 1) * sizeof(int));
    
    for (int i = 1; i <= n; i++) {
        int is_sink = 1;
        for (int j = 1; j <= n; j++) {
            if (graph[i][j] == 1) {
                is_sink = 0;
                break;
            }
        }
        if (is_sink) {
            sinks[count++] = i;
        }
    }
    
    // Output the number of sinks
    printf("%d\n", count);
    
    // Output the sink nodes
    for (int i = 0; i < count; i++) {
        printf("%d ", sinks[i]);
    }
    printf("\n");
    
    // Free memory
    for (int i = 0; i <= n; i++) {
        free(graph[i]);
    }
    free(graph);
    free(sinks);
    
    return 0;
}
```

## Explanation

1. **Input Reading**: Read the number of nodes `n` and edges `m`
2. **Graph Representation**: Use a 2D boolean matrix to represent the adjacency list
3. **Edge Population**: For each edge (u,v), mark `graph[u][v] = 1`
4. **Sink Detection**: For each node i, check if all entries in row i are 0 (no outgoing edges)
5. **Output**: Print count of sinks followed by the sink nodes

## Time Complexity
- O(n²) for checking all nodes and their outgoing edges
- Space complexity: O(n²) for the adjacency matrix

## Sample Input/Output

**Input:**
```
4 4
1 2
2 3
3 4
4 1
```

**Output:**
```
0
```

**Input:**
```
4 3
1 2
2 3
3 4
```

**Output:**
```
1
4
```