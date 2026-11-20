# Rosalind Problem: Shortest Cycle Through a Given Edge

## Problem Understanding

Given an undirected graph and a specific edge, we need to find the shortest cycle that includes this particular edge.

## Approach

1. Remove the given edge from the graph
2. Find the shortest path between the two vertices of the removed edge in the remaining graph
3. Add the removed edge back to form the cycle
4. The shortest cycle length = shortest path + 1 (the removed edge)

## Solution Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define MAX_NODES 1000
#define INF INT_MAX

// Graph representation using adjacency list
typedef struct Node {
    int vertex;
    struct Node* next;
} Node;

typedef struct {
    Node* head[MAX_NODES];
    int num_vertices;
    int num_edges;
} Graph;

// Queue for BFS
typedef struct {
    int items[MAX_NODES];
    int front;
    int rear;
} Queue;

// Function to create a new node
Node* create_node(int vertex) {
    Node* new_node = (Node*)malloc(sizeof(Node));
    new_node->vertex = vertex;
    new_node->next = NULL;
    return new_node;
}

// Function to create a graph
Graph* create_graph(int vertices) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->num_vertices = vertices;
    graph->num_edges = 0;
    
    for (int i = 0; i < MAX_NODES; i++) {
        graph->head[i] = NULL;
    }
    
    return graph;
}

// Function to add edge to graph
void add_edge(Graph* graph, int src, int dest) {
    // Add edge from src to dest
    Node* new_node = create_node(dest);
    new_node->next = graph->head[src];
    graph->head[src] = new_node;
    
    // Add edge from dest to src (undirected)
    new_node = create_node(src);
    new_node->next = graph->head[dest];
    graph->head[dest] = new_node;
    
    graph->num_edges++;
}

// Function to initialize queue
Queue* create_queue() {
    Queue* queue = (Queue*)malloc(sizeof(Queue));
    queue->front = -1;
    queue->rear = -1;
    return queue;
}

// Function to check if queue is empty
int is_empty(Queue* queue) {
    return queue->front == -1;
}

// Function to enqueue
void enqueue(Queue* queue, int item) {
    if (queue->rear == MAX_NODES - 1) {
        printf("Queue overflow\n");
        return;
    }
    
    if (queue->front == -1) {
        queue->front = 0;
    }
    
    queue->rear++;
    queue->items[queue->rear] = item;
}

// Function to dequeue
int dequeue(Queue* queue) {
    if (is_empty(queue)) {
        printf("Queue underflow\n");
        return -1;
    }
    
    int item = queue->items[queue->front];
    
    if (queue->front == queue->rear) {
        queue->front = -1;
        queue->rear = -1;
    } else {
        queue->front++;
    }
    
    return item;
}

// Function to find shortest path between two vertices using BFS
int bfs_shortest_path(Graph* graph, int start, int end) {
    if (start == end) return 0;
    
    int visited[MAX_NODES] = {0};
    int distance[MAX_NODES] = {0};
    
    Queue* queue = create_queue();
    
    visited[start] = 1;
    distance[start] = 0;
    enqueue(queue, start);
    
    while (!is_empty(queue)) {
        int current = dequeue(queue);
        
        Node* temp = graph->head[current];
        while (temp != NULL) {
            int neighbor = temp->vertex;
            
            if (!visited[neighbor]) {
                visited[neighbor] = 1;
                distance[neighbor] = distance[current] + 1;
                
                if (neighbor == end) {
                    return distance[neighbor];
                }
                
                enqueue(queue, neighbor);
            }
            temp = temp->next;
        }
    }
    
    return -1; // No path found
}

// Function to find shortest cycle through a given edge
int shortest_cycle_through_edge(Graph* graph, int u, int v) {
    // Remove the edge (u, v)
    Node* prev = NULL;
    Node* current = graph->head[u];
    
    while (current != NULL) {
        if (current->vertex == v) {
            if (prev == NULL) {
                graph->head[u] = current->next;
            } else {
                prev->next = current->next;
            }
            free(current);
            break;
        }
        prev = current;
        current = current->next;
    }
    
    // Remove the corresponding edge from v to u
    prev = NULL;
    current = graph->head[v];
    
    while (current != NULL) {
        if (current->vertex == u) {
            if (prev == NULL) {
                graph->head[v] = current->next;
            } else {
                prev->next = current->next;
            }
            free(current);
            break;
        }
        prev = current;
        current = current->next;
    }
    
    // Find shortest path between u and v in the modified graph
    int shortest_path_length = bfs_shortest_path(graph, u, v);
    
    // Restore the edge (u, v)
    add_edge(graph, u, v);
    
    // Return cycle length (shortest path + 1 for the removed edge)
    if (shortest_path_length == -1) {
        return -1; // No cycle possible
    }
    
    return shortest_path_length + 1;
}

// Function to print graph (for debugging)
void print_graph(Graph* graph) {
    printf("Graph:\n");
    for (int i = 0; i < graph->num_vertices; i++) {
        printf("%d: ", i);
        Node* temp = graph->head[i];
        while (temp != NULL) {
            printf("%d ", temp->vertex);
            temp = temp->next;
        }
        printf("\n");
    }
}

int main() {
    // Read input
    int n, m;
    scanf("%d %d", &n, &m);
    
    Graph* graph = create_graph(n);
    
    // Read edges
    for (int i = 0; i < m; i++) {
        int u, v;
        scanf("%d %d", &u, &v);
        add_edge(graph, u, v);
    }
    
    // Read the specific edge
    int u, v;
    scanf("%d %d", &u, &v);
    
    // Find shortest cycle through given edge
    int result = shortest_cycle_through_edge(graph, u, v);
    
    if (result == -1) {
        printf("No cycle exists\n");
    } else {
        printf("%d\n", result);
    }
    
    return 0;
}
```

## Explanation

The solution works as follows:

1. **Graph Representation**: Uses adjacency list representation for efficient storage and traversal
2. **Edge Removal**: Temporarily removes the specified edge from the graph
3. **BFS Search**: Uses breadth-first search to find the shortest path between the two vertices of the removed edge
4. **Cycle Calculation**: The length of the shortest cycle is the shortest path length plus 1 (the removed edge)
5. **Edge Restoration**: Restores the removed edge for potential future operations

## Time Complexity
- O(V + E) where V is the number of vertices and E is the number of edges
- BFS traversal takes O(V + E) time
- Edge removal and restoration are O(1) operations

## Space Complexity
- O(V + E) for storing the graph and auxiliary data structures

The algorithm correctly handles the case where no cycle exists through the given edge, returning -1 in such cases.