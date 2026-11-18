# Rosalind Problem: Testing Acyclicity - C Solution

## Problem Understanding

The problem asks us to determine whether a directed graph is acyclic (i.e., contains no cycles). This is equivalent to checking if we can perform a topological sort on the graph.

## Approach

We'll use **Kahn's algorithm** for topological sorting:
1. Calculate in-degrees of all vertices
2. Add all vertices with in-degree 0 to a queue
3. Process vertices from queue, reducing in-degrees of neighbors
4. If we process all vertices, the graph is acyclic; otherwise, it contains cycles

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NODES 10000

// Structure for adjacency list
typedef struct Node {
    int data;
    struct Node* next;
} Node;

// Structure for adjacency list
typedef struct {
    Node* head;
} AdjList;

// Function to create a new node
Node* createNode(int data) {
    Node* newNode = (Node*)malloc(sizeof(Node));
    newNode->data = data;
    newNode->next = NULL;
    return newNode;
}

// Function to add edge to adjacency list
void addEdge(AdjList* adjList, int src, int dest) {
    Node* newNode = createNode(dest);
    newNode->next = adjList[src].head;
    adjList[src].head = newNode;
}

// Function to check if graph is acyclic
int isAcyclic(int n, AdjList* adjList) {
    // Calculate in-degrees
    int* inDegree = (int*)calloc(n + 1, sizeof(int));
    
    for (int i = 1; i <= n; i++) {
        Node* current = adjList[i].head;
        while (current != NULL) {
            inDegree[current->data]++;
            current = current->next;
        }
    }
    
    // Initialize queue with vertices having in-degree 0
    int* queue = (int*)malloc((n + 1) * sizeof(int));
    int front = 0, rear = 0;
    
    for (int i = 1; i <= n; i++) {
        if (inDegree[i] == 0) {
            queue[rear++] = i;
        }
    }
    
    int processedCount = 0;
    
    // Process vertices in topological order
    while (front < rear) {
        int current = queue[front++];
        processedCount++;
        
        // Reduce in-degree of neighbors
        Node* currentEdge = adjList[current].head;
        while (currentEdge != NULL) {
            inDegree[currentEdge->data]--;
            if (inDegree[currentEdge->data] == 0) {
                queue[rear++] = currentEdge->data;
            }
            currentEdge = currentEdge->next;
        }
    }
    
    free(inDegree);
    free(queue);
    
    // If we processed all vertices, graph is acyclic
    return (processedCount == n) ? 0 : 1;
}

int main() {
    int n, m;
    scanf("%d %d", &n, &m);
    
    // Create adjacency list
    AdjList* adjList = (AdjList*)malloc((n + 1) * sizeof(AdjList));
    
    // Initialize adjacency lists
    for (int i = 0; i <= n; i++) {
        adjList[i].head = NULL;
    }
    
    // Read edges
    int src, dest;
    for (int i = 0; i < m; i++) {
        scanf("%d %d", &src, &dest);
        addEdge(adjList, src, dest);
    }
    
    // Check if graph is acyclic
    int result = isAcyclic(n, adjList);
    
    printf("%d\n", result);
    
    // Free memory
    for (int i = 1; i <= n; i++) {
        Node* current = adjList[i].head;
        while (current != NULL) {
            Node* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(adjList);
    
    return 0;
}
```

## Explanation

1. **Data Structures**: 
   - We use an adjacency list representation with linked lists
   - In-degree array to track incoming edges for each vertex

2. **Algorithm Steps**:
   - Calculate in-degrees of all vertices by traversing the adjacency list
   - Initialize queue with all vertices having in-degree 0
   - Process vertices from queue, reducing in-degrees of neighbors
   - If we can process all vertices (topological sort), graph is acyclic

3. **Time Complexity**: O(V + E) where V is vertices and E is edges
4. **Space Complexity**: O(V + E) for adjacency list and auxiliary arrays

## Sample Input/Output

**Input:**
```
4 4
1 2
2 3
3 1
4 1
```

**Output:**
```
1
```

The output `1` indicates the graph contains a cycle (1→2→3→1), making it cyclic.

**Input:**
```
4 3
1 2
2 3
3 4
```

**Output:**
```
0
```

The output `0` indicates the graph is acyclic (topological sort possible).