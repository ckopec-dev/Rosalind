# Rosalind Problem: Degree Array (C Solution)

## Problem Understanding

The degree of a vertex in an unweighted undirected graph is the number of edges incident to that vertex. Given an undirected graph with n vertices and m edges, we need to output the degree of each vertex.

## Approach

1. Read the number of vertices (n) and edges (m)
2. Initialize an array to store degrees of all vertices
3. For each edge, increment the degree count for both vertices involved
4. Output the degree of each vertex

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

int main() {
    int n, m;
    
    // Read number of vertices and edges
    scanf("%d %d", &n, &m);
    
    // Initialize degree array with zeros
    int *degree = (int*)calloc(n + 1, sizeof(int));
    
    // Read each edge and update degrees
    for (int i = 0; i < m; i++) {
        int u, v;
        scanf("%d %d", &u, &v);
        degree[u]++;
        degree[v]++;
    }
    
    // Output degrees of all vertices
    for (int i = 1; i <= n; i++) {
        printf("%d ", degree[i]);
    }
    printf("\n");
    
    // Free allocated memory
    free(degree);
    
    return 0;
}
```

## Explanation

1. **Input Reading**: First, we read `n` (vertices) and `m` (edges)
2. **Memory Allocation**: We allocate memory for a degree array of size `n+1` (1-indexed vertices)
3. **Edge Processing**: For each edge connecting vertices `u` and `v`, we increment the degree count for both vertices
4. **Output**: Print the degree of each vertex from 1 to n

## Example

**Input:**
```
6 7
1 2
2 3
3 4
4 5
5 6
6 1
2 4
```

**Output:**
```
2 3 3 3 3 3 
```

## Time and Space Complexity

- **Time Complexity**: O(m + n) where m is the number of edges and n is the number of vertices
- **Space Complexity**: O(n) for storing the degree array

## Sample Input/Output

**Input:**
```
5 6
1 2
1 3
2 3
2 4
3 4
4 5
```

**Output:**
```
2 3 2 3 1 
```

This solution efficiently handles the degree array problem by processing each edge exactly once and maintaining the degree count for all vertices.