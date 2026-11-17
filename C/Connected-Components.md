# Rosalind Problem: Connected Components in C

## Problem Understanding

The connected components problem asks us to find how many connected components exist in an undirected graph. A connected component is a subgraph where all vertices are connected to each other through paths, and which is disconnected from the rest of the graph.

## Solution Approach

I'll use a Union-Find (Disjoint Set) data structure with path compression and union by rank to efficiently solve this problem:

1. Initialize each vertex as its own parent
2. Process all edges and perform union operations
3. Count the number of distinct root nodes (connected components)

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

// Union-Find data structure
typedef struct {
    int *parent;
    int *rank;
    int n;
} UnionFind;

// Initialize Union-Find structure
UnionFind* init_union_find(int n) {
    UnionFind *uf = (UnionFind*)malloc(sizeof(UnionFind));
    uf->n = n;
    uf->parent = (int*)malloc(n * sizeof(int));
    uf->rank = (int*)malloc(n * sizeof(int));
    
    // Initialize each node as its own parent
    for (int i = 0; i < n; i++) {
        uf->parent[i] = i;
        uf->rank[i] = 0;
    }
    
    return uf;
}

// Find root with path compression
int find(UnionFind *uf, int x) {
    if (uf->parent[x] != x) {
        uf->parent[x] = find(uf, uf->parent[x]); // Path compression
    }
    return uf->parent[x];
}

// Union two sets
void union_sets(UnionFind *uf, int x, int y) {
    int root_x = find(uf, x);
    int root_y = find(uf, y);
    
    if (root_x != root_y) {
        // Union by rank
        if (uf->rank[root_x] < uf->rank[root_y]) {
            uf->parent[root_x] = root_y;
        } else if (uf->rank[root_x] > uf->rank[root_y]) {
            uf->parent[root_y] = root_x;
        } else {
            uf->parent[root_y] = root_x;
            uf->rank[root_x]++;
        }
    }
}

// Count connected components
int count_components(UnionFind *uf) {
    int count = 0;
    for (int i = 0; i < uf->n; i++) {
        if (uf->parent[i] == i) {
            count++;
        }
    }
    return count;
}

// Free memory
void free_union_find(UnionFind *uf) {
    free(uf->parent);
    free(uf->rank);
    free(uf);
}

int main() {
    int n, m;
    
    // Read number of vertices and edges
    scanf("%d %d", &n, &m);
    
    // Initialize Union-Find structure
    UnionFind *uf = init_union_find(n);
    
    // Process all edges
    for (int i = 0; i < m; i++) {
        int u, v;
        scanf("%d %d", &u, &v);
        // Convert to 0-based indexing
        union_sets(uf, u - 1, v - 1);
    }
    
    // Count connected components
    int components = count_components(uf);
    printf("%d\n", components);
    
    // Free memory
    free_union_find(uf);
    
    return 0;
}
```

## How It Works

1. **Union-Find Structure**: 
   - `parent[]`: Tracks parent of each node
   - `rank[]`: Tracks depth/size of tree for union optimization

2. **Find Operation**: 
   - Uses path compression to make future searches faster
   - Returns the root of the set containing element x

3. **Union Operation**: 
   - Uses union by rank to keep trees balanced
   - Merges two sets if they're different

4. **Component Counting**: 
   - Counts how many nodes are their own parent (roots)
   - Each root represents one connected component

## Input Format
```
n m
u1 v1
u2 v2
...
um vm
```

Where:
- n = number of vertices
- m = number of edges
- Each line after the first contains an edge between vertices ui and vi

## Time Complexity
- **Time**: O(m × α(n)) where α is the inverse Ackermann function (practically constant)
- **Space**: O(n)

## Example

Input:
```
4 2
1 2
3 4
```

Output:
```
2
```

This represents a graph with 4 vertices and 2 edges, forming 2 connected components: {1,2} and {3,4}.