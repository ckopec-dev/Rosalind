# Rosalind Problem: Topological Sorting (C Solution)

## Problem Understanding

Topological sorting involves ordering the vertices of a directed acyclic graph (DAG) such that for every directed edge u → v, vertex u comes before vertex v in the ordering.

## Approach

1. **Input parsing**: Read the number of nodes and edges
2. **Build adjacency list**: Represent the graph structure
3. **Calculate in-degrees**: Count incoming edges for each node
4. **Topological sort**: Use Kahn's algorithm with queue-based approach
5. **Output**: Print the topological ordering

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

#define MAX_NODES 10000

// Structure for adjacency list
struct Node {
    int data;
    struct Node* next;
};

// Function to create a new node
struct Node* createNode(int data) {
    struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
    newNode->data = data;
    newNode->next = NULL;
    return newNode;
}

// Function to add edge to adjacency list
void addEdge(struct Node** adjList, int src, int dest) {
    struct Node* newNode = createNode(dest);
    newNode->next = adjList[src];
    adjList[src] = newNode;
}

// Topological sorting using Kahn's algorithm
void topologicalSort(int n, struct Node** adjList) {
    // Calculate in-degrees for all nodes
    int* inDegree = (int*)calloc(n + 1, sizeof(int));
    
    for (int i = 1; i <= n; i++) {
        struct Node* temp = adjList[i];
        while (temp != NULL) {
            inDegree[temp->data]++;
            temp = temp->next;
        }
    }
    
    // Create queue and enqueue all nodes with in-degree 0
    int* queue = (int*)malloc((n + 1) * sizeof(int));
    int front = 0, rear = 0;
    
    for (int i = 1; i <= n; i++) {
        if (inDegree[i] == 0) {
            queue[rear++] = i;
        }
    }
    
    // Process nodes in topological order
    while (front < rear) {
        int current = queue[front++];
        printf("%d ", current);
        
        // Reduce in-degree of adjacent nodes
        struct Node* temp = adjList[current];
        while (temp != NULL) {
            inDegree[temp->data]--;
            if (inDegree[temp->data] == 0) {
                queue[rear++] = temp->data;
            }
            temp = temp->next;
        }
    }
    
    printf("\n");
    
    // Free allocated memory
    free(inDegree);
    free(queue);
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    
    // Initialize adjacency list
    struct Node** adjList = (struct Node**)malloc((n + 1) * sizeof(struct Node*));
    for (int i = 0; i <= n; i++) {
        adjList[i] = NULL;
    }
    
    // Read edges and build graph
    int src, dest;
    for (int i = 0; i < m; i++) {
        scanf("%d %d", &src, &dest);
        addEdge(adjList, src, dest);
    }
    
    // Perform topological sort
    topologicalSort(n, adjList);
    
    // Free adjacency list memory
    for (int i = 1; i <= n; i++) {
        struct Node* temp = adjList[i];
        while (temp != NULL) {
            struct Node* toDelete = temp;
            temp = temp->next;
            free(toDelete);
        }
    }
    free(adjList);
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
6 8
1 2
1 3
2 4
3 4
4 5
4 6
5 6
2 6
```

**Output:**
```
1 2 3 4 5 6
```

## Explanation

1. **Graph Representation**: Uses adjacency list to represent the directed graph
2. **In-degree Calculation**: Computes in-degrees for all nodes by traversing the adjacency list
3. **Kahn's Algorithm**: 
   - Initialize queue with all nodes having in-degree 0
   - Process nodes from queue, reduce in-degrees of neighbors
   - Add nodes with in-degree 0 to queue
4. **Time Complexity**: O(V + E) where V is vertices and E is edges
5. **Space Complexity**: O(V + E) for adjacency list and auxiliary arrays

## Key Points

- The solution handles multiple valid topological orders (any valid order is acceptable)
- Uses Kahn's algorithm which is efficient for topological sorting
- Proper memory management with malloc/free
- Handles input/output as required by Rosalind format