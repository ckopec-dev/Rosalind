# Rosalind Problem: Semi-Connected Graph

## Problem Understanding

A directed graph is semi-connected if for every pair of vertices u and v, there is either a path from u to v or a path from v to u (or both).

## Approach

To solve this problem:
1. Find the strongly connected components (SCCs) using Kosaraju's algorithm
2. Create a condensation graph (component graph) where each SCC is a node
3. Check if the condensation graph forms a directed path (i.e., it's a DAG with a unique topological ordering)
4. If yes, the original graph is semi-connected

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 10005

// Global variables
int adj[MAXN][MAXN];
int adj_rev[MAXN][MAXN];
int visited[MAXN];
int scc[MAXN];
int stack[MAXN];
int stack_size;
int scc_count;
int in_degree[MAXN];
int out_degree[MAXN];

// DFS for first pass (forward graph)
void dfs1(int u) {
    visited[u] = 1;
    for (int v = 1; v <= MAXN; v++) {
        if (adj[u][v] && !visited[v]) {
            dfs1(v);
        }
    }
    stack[++stack_size] = u;
}

// DFS for second pass (reverse graph)
void dfs2(int u) {
    visited[u] = 1;
    scc[u] = scc_count;
    for (int v = 1; v <= MAXN; v++) {
        if (adj_rev[u][v] && !visited[v]) {
            dfs2(v);
        }
    }
}

// Kosaraju's algorithm to find SCCs
void kosaraju(int n) {
    // Reset visited array
    memset(visited, 0, sizeof(visited));
    
    // First pass - get finishing times
    stack_size = 0;
    for (int i = 1; i <= n; i++) {
        if (!visited[i]) {
            dfs1(i);
        }
    }
    
    // Second pass - find SCCs
    memset(visited, 0, sizeof(visited));
    scc_count = 0;
    
    while (stack_size > 0) {
        int u = stack[stack_size--];
        if (!visited[u]) {
            scc_count++;
            dfs2(u);
        }
    }
}

// Check if semi-connected
int is_semi_connected(int n) {
    // Reset degrees
    memset(in_degree, 0, sizeof(in_degree));
    memset(out_degree, 0, sizeof(out_degree));
    
    // Count in and out degrees of SCCs
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (adj[i][j] && scc[i] != scc[j]) {
                in_degree[scc[j]]++;
                out_degree[scc[i]]++;
            }
        }
    }
    
    // Count SCCs with in_degree = 0 and out_degree = 0
    int zero_in = 0, zero_out = 0;
    for (int i = 1; i <= scc_count; i++) {
        if (in_degree[i] == 0) zero_in++;
        if (out_degree[i] == 0) zero_out++;
    }
    
    // For semi-connected graph, there should be exactly one node with in_degree=0 
    // and exactly one node with out_degree=0, and all other nodes have both in and out degree > 0
    // This means we can form a directed path through all SCCs
    
    int total_nodes = scc_count;
    
    // If there's only one SCC, it's semi-connected
    if (total_nodes == 1) return 1;
    
    // Check if condensation graph is a directed path
    int start_count = 0, end_count = 0;
    for (int i = 1; i <= scc_count; i++) {
        if (in_degree[i] == 0) start_count++;
        if (out_degree[i] == 0) end_count++;
    }
    
    // For a path, there should be exactly one start and one end
    if (start_count == 1 && end_count == 1) {
        // Check that we have a valid path (in_degree + out_degree = total_nodes - 1 for each intermediate node)
        int valid_path = 1;
        for (int i = 1; i <= scc_count; i++) {
            if (in_degree[i] == 0 && out_degree[i] == 0) continue; // This is handled in the above check
            if (i != 1 && i != total_nodes && (in_degree[i] + out_degree[i] != 2)) {
                valid_path = 0;
                break;
            }
        }
        return valid_path;
    }
    
    return 0;
}

// Simplified approach for semi-connected check
int is_semi_connected_simple(int n) {
    // Run Kosaraju's algorithm to find SCCs
    kosaraju(n);
    
    // Create component graph
    int component_adj[MAXN][MAXN];
    memset(component_adj, 0, sizeof(component_adj));
    
    // For each edge in original graph, if it connects different SCCs, 
    // add an edge in the component graph
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (adj[i][j] && scc[i] != scc[j]) {
                component_adj[scc[i]][scc[j]] = 1;
            }
        }
    }
    
    // Check if component graph forms a directed path
    int in_degree[MAXN] = {0};
    int out_degree[MAXN] = {0};
    
    for (int i = 1; i <= scc_count; i++) {
        for (int j = 1; j <= scc_count; j++) {
            if (component_adj[i][j]) {
                in_degree[j]++;
                out_degree[i]++;
            }
        }
    }
    
    // Count nodes with in/out degree 0
    int start_nodes = 0, end_nodes = 0;
    for (int i = 1; i <= scc_count; i++) {
        if (in_degree[i] == 0) start_nodes++;
        if (out_degree[i] == 0) end_nodes++;
    }
    
    // For a directed path: exactly one start and one end
    if (start_nodes == 1 && end_nodes == 1) {
        // Check that it's a valid path
        int path_length = scc_count;
        for (int i = 1; i <= scc_count; i++) {
            if (in_degree[i] + out_degree[i] != 2 && 
                (i != 1 && i != scc_count)) {
                // Intermediate nodes should have in_degree + out_degree = 2
                return 0;
            }
        }
        return 1;
    }
    
    return 0;
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    
    // Reset adjacency matrices
    memset(adj, 0, sizeof(adj));
    memset(adj_rev, 0, sizeof(adj_rev));
    
    for (int i = 0; i < m; i++) {
        int u, v;
        scanf("%d %d", &u, &v);
        adj[u][v] = 1;
        adj_rev[v][u] = 1;
    }
    
    // Check if the graph is semi-connected
    kosaraju(n);
    
    // If there's only one SCC, it's semi-connected
    if (scc_count == 1) {
        printf("1\n");
        return 0;
    }
    
    // Build component graph
    int comp_adj[MAXN][MAXN];
    memset(comp_adj, 0, sizeof(comp_adj));
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (adj[i][j] && scc[i] != scc[j]) {
                comp_adj[scc[i]][scc[j]] = 1;
            }
        }
    }
    
    // Calculate in and out degrees of components
    int comp_in[MAXN] = {0};
    int comp_out[MAXN] = {0};
    
    for (int i = 1; i <= scc_count; i++) {
        for (int j = 1; j <= scc_count; j++) {
            if (comp_adj[i][j]) {
                comp_in[j]++;
                comp_out[i]++;
            }
        }
    }
    
    // Check if it's a path
    int start_count = 0, end_count = 0;
    for (int i = 1; i <= scc_count; i++) {
        if (comp_in[i] == 0) start_count++;
        if (comp_out[i] == 0) end_count++;
    }
    
    // For a semi-connected graph, we need exactly one start and one end
    int is_semi = 0;
    if (start_count == 1 && end_count == 1) {
        is_semi = 1;
        // Validate that it's actually a path
        for (int i = 1; i <= scc_count; i++) {
            if (i != 1 && i != scc_count && (comp_in[i] + comp_out[i] != 2)) {
                is_semi = 0;
                break;
            }
        }
    }
    
    printf("%d\n", is_semi);
    
    return 0;
}
```

## Explanation

The solution uses the following approach:

1. **Kosaraju's Algorithm**: Finds all strongly connected components (SCCs) in the directed graph
2. **Component Graph Construction**: Creates a condensed graph where each node represents an SCC
3. **Path Validation**: Checks if the component graph forms a directed path:
   - Exactly one node with in-degree 0 (start)
   - Exactly one node with out-degree 0 (end)  
   - All intermediate nodes have in-degree + out-degree = 2

## Time Complexity
- O(V + E) for Kosaraju's algorithm
- O(V²) for building component graph
- Overall: O(V² + E)

## Space Complexity
- O(V²) for adjacency matrices
- O(V) for auxiliary arrays

This solution correctly identifies whether a directed graph is semi-connected by checking if its condensation graph forms a directed path.