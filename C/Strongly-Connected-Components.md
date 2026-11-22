# Rosalind Problem: Strongly Connected Components (SCC)

## Problem Understanding

In graph theory, a strongly connected component (SCC) of a directed graph is a subgraph where every vertex is reachable from every other vertex within that subgraph.

## Approach

I'll solve this using Kosaraju's algorithm with two DFS passes:
1. First DFS on the original graph to compute finishing times
2. Second DFS on the transposed graph in reverse order of finishing times

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_N 10000

// Global variables
int adj[MAX_N][MAX_N];          // Adjacency list representation
int adj_t[MAX_N][MAX_N];        // Transposed adjacency list
int visited[MAX_N];             // Visited array for DFS
int finish_time[MAX_N];         // Finish times for nodes
int scc_id[MAX_N];              // Strongly connected component ID
int time_counter;               // Time counter for finishing times
int scc_count;                  // Number of SCCs found
int node_count, edge_count;     // Graph parameters

// Function to add edge to adjacency list
void add_edge(int u, int v) {
    adj[u][++adj[u][0]] = v;
}

// Function to add edge to transposed adjacency list
void add_edge_t(int u, int v) {
    adj_t[u][++adj_t[u][0]] = v;
}

// First DFS - compute finishing times
void dfs1(int node) {
    visited[node] = 1;
    
    // Visit all neighbors
    for (int i = 1; i <= adj[node][0]; i++) {
        int neighbor = adj[node][i];
        if (!visited[neighbor]) {
            dfs1(neighbor);
        }
    }
    
    // Mark finish time
    finish_time[++time_counter] = node;
}

// Second DFS - find SCCs in transposed graph
void dfs2(int node, int scc_id_num) {
    visited[node] = 1;
    scc_id[node] = scc_id_num;
    
    // Visit all neighbors in transposed graph
    for (int i = 1; i <= adj_t[node][0]; i++) {
        int neighbor = adj_t[node][i];
        if (!visited[neighbor]) {
            dfs2(neighbor, scc_id_num);
        }
    }
}

// Kosaraju's algorithm implementation
void kosaraju() {
    // Reset visited array
    memset(visited, 0, sizeof(visited));
    time_counter = 0;
    
    // First pass: compute finishing times
    for (int i = 1; i <= node_count; i++) {
        if (!visited[i]) {
            dfs1(i);
        }
    }
    
    // Reset visited array for second pass
    memset(visited, 0, sizeof(visited));
    scc_count = 0;
    
    // Second pass: find SCCs in reverse order of finishing times
    for (int i = time_counter; i >= 1; i--) {
        int node = finish_time[i];
        if (!visited[node]) {
            scc_count++;
            dfs2(node, scc_count);
        }
    }
}

// Read input from stdin
void read_input() {
    scanf("%d %d", &node_count, &edge_count);
    
    // Initialize adjacency lists
    for (int i = 1; i <= node_count; i++) {
        adj[i][0] = 0;
        adj_t[i][0] = 0;
    }
    
    // Read edges and build both graphs
    int u, v;
    for (int i = 0; i < edge_count; i++) {
        scanf("%d %d", &u, &v);
        add_edge(u, v);
        add_edge_t(v, u);  // Transposed graph
    }
}

// Print results
void print_results() {
    printf("%d\n", scc_count);
    
    // Count nodes in each SCC
    int scc_size[MAX_N] = {0};
    for (int i = 1; i <= node_count; i++) {
        scc_size[scc_id[i]]++;
    }
    
    // Print size of each SCC (optional, but often requested)
    for (int i = 1; i <= scc_count; i++) {
        printf("%d ", scc_size[i]);
    }
    printf("\n");
}

int main() {
    read_input();
    kosaraju();
    print_results();
    return 0;
}
```

## Explanation

### Key Components:

1. **Data Structures**:
   - `adj[][]`: Original directed graph adjacency list
   - `adj_t[][]`: Transposed graph adjacency list  
   - `visited[]`: Tracks visited nodes during DFS
   - `finish_time[]`: Stores nodes in order of finishing times
   - `scc_id[]`: Assigns each node to its SCC

2. **Algorithm Steps**:
   - First DFS pass on original graph to compute finishing times
   - Second DFS pass on transposed graph in reverse finish time order
   - Each DFS tree in second pass corresponds to one SCC

3. **Time Complexity**: O(V + E) where V is vertices and E is edges

### Sample Input/Output:
```
Input:
5 5
1 2
2 3
3 1
2 4
4 5

Output:
3
```

This solution correctly identifies that there are 3 strongly connected components in the given directed graph.