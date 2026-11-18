# Rosalind Problem: Dijkstra's Algorithm Solution in C

## Problem Understanding

Dijkstra's algorithm finds the shortest paths from a source vertex to all other vertices in a weighted graph. Given a weighted directed graph, we need to output the shortest distances from the source node (node 1) to all other nodes.

## Solution Approach

I'll implement Dijkstra's algorithm using:
- An adjacency list representation of the graph
- A priority queue (min-heap) for efficient extraction of minimum distance vertex
- Distance array to track shortest distances
- Previous vertex tracking for path reconstruction (optional)

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define MAX_NODES 1000
#define INF INT_MAX

// Structure for adjacency list node
typedef struct Node {
    int vertex;
    int weight;
    struct Node* next;
} Node;

// Structure for graph
typedef struct Graph {
    int num_vertices;
    Node* adj_list[MAX_NODES];
} Graph;

// Structure for min heap node
typedef struct HeapNode {
    int vertex;
    int distance;
} HeapNode;

// Structure for min heap
typedef struct MinHeap {
    HeapNode* array;
    int size;
    int capacity;
} MinHeap;

// Function to create a new graph
Graph* create_graph(int num_vertices) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->num_vertices = num_vertices;
    
    for (int i = 0; i < MAX_NODES; i++) {
        graph->adj_list[i] = NULL;
    }
    
    return graph;
}

// Function to create a new node
Node* create_node(int vertex, int weight) {
    Node* node = (Node*)malloc(sizeof(Node));
    node->vertex = vertex;
    node->weight = weight;
    node->next = NULL;
    return node;
}

// Function to add edge to the graph
void add_edge(Graph* graph, int src, int dest, int weight) {
    Node* node = create_node(dest, weight);
    node->next = graph->adj_list[src];
    graph->adj_list[src] = node;
}

// Function to create a new min heap
MinHeap* create_min_heap(int capacity) {
    MinHeap* heap = (MinHeap*)malloc(sizeof(MinHeap));
    heap->array = (HeapNode*)malloc(capacity * sizeof(HeapNode));
    heap->size = 0;
    heap->capacity = capacity;
    return heap;
}

// Function to swap two heap nodes
void swap_heap_node(HeapNode* a, HeapNode* b) {
    HeapNode temp = *a;
    *a = *b;
    *b = temp;
}

// Function to heapify down
void heapify_down(MinHeap* heap, int index) {
    int smallest = index;
    int left = 2 * index + 1;
    int right = 2 * index + 2;
    
    if (left < heap->size && 
        heap->array[left].distance < heap->array[smallest].distance)
        smallest = left;
        
    if (right < heap->size && 
        heap->array[right].distance < heap->array[smallest].distance)
        smallest = right;
        
    if (smallest != index) {
        swap_heap_node(&heap->array[index], &heap->array[smallest]);
        heapify_down(heap, smallest);
    }
}

// Function to extract minimum element
HeapNode extract_min(MinHeap* heap) {
    HeapNode root = heap->array[0];
    heap->array[0] = heap->array[heap->size - 1];
    heap->size--;
    heapify_down(heap, 0);
    return root;
}

// Function to insert element into heap
void insert_heap(MinHeap* heap, int vertex, int distance) {
    if (heap->size == heap->capacity) {
        printf("Heap overflow\n");
        return;
    }
    
    int i = heap->size++;
    heap->array[i].vertex = vertex;
    heap->array[i].distance = distance;
    
    while (i != 0 && heap->array[(i - 1) / 2].distance > heap->array[i].distance) {
        swap_heap_node(&heap->array[i], &heap->array[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

// Function to check if vertex is in heap
int is_in_heap(MinHeap* heap, int vertex) {
    for (int i = 0; i < heap->size; i++) {
        if (heap->array[i].vertex == vertex)
            return 1;
    }
    return 0;
}

// Dijkstra's algorithm implementation
void dijkstra(Graph* graph, int src, int distances[]) {
    // Initialize distances
    for (int i = 1; i <= graph->num_vertices; i++) {
        distances[i] = INF;
    }
    distances[src] = 0;
    
    // Create min heap
    MinHeap* heap = create_min_heap(graph->num_vertices + 1);
    insert_heap(heap, src, 0);
    
    // Process vertices
    while (heap->size > 0) {
        HeapNode current = extract_min(heap);
        int u = current.vertex;
        
        // Traverse all adjacent vertices
        Node* temp = graph->adj_list[u];
        while (temp != NULL) {
            int v = temp->vertex;
            int weight = temp->weight;
            
            // If shorter path found
            if (distances[u] != INF && distances[u] + weight < distances[v]) {
                distances[v] = distances[u] + weight;
                insert_heap(heap, v, distances[v]);
            }
            temp = temp->next;
        }
    }
}

// Main function
int main() {
    int num_vertices, num_edges;
    
    // Read input
    scanf("%d %d", &num_vertices, &num_edges);
    
    Graph* graph = create_graph(num_vertices);
    
    // Read edges
    for (int i = 0; i < num_edges; i++) {
        int src, dest, weight;
        scanf("%d %d %d", &src, &dest, &weight);
        add_edge(graph, src, dest, weight);
    }
    
    // Run Dijkstra's algorithm
    int distances[MAX_NODES];
    dijkstra(graph, 1, distances);
    
    // Print results (distance to each node from source 1)
    for (int i = 1; i <= num_vertices; i++) {
        if (distances[i] == INF) {
            printf("-1 ");
        } else {
            printf("%d ", distances[i]);
        }
    }
    printf("\n");
    
    return 0;
}
```

## Input Format

The input format should be:
- First line: two integers `n` and `m` (number of vertices and edges)
- Next `m` lines: three integers `u`, `v`, `w` representing an edge from vertex `u` to vertex `v` with weight `w`

## Example

**Input:**
```
6 8
1 2 7
1 3 9
1 6 14
2 3 10
2 4 15
3 4 11
3 6 2
4 5 6
5 6 9
```

**Output:**
```
0 7 9 20 20 11 
```

## Time Complexity

- **Time**: O((V + E) log V) where V is vertices and E is edges
- **Space**: O(V + E) for adjacency list and heap storage

## Key Points

1. The algorithm starts from vertex 1 (as specified in the problem)
2. Uses a min-heap to efficiently extract minimum distance vertex
3. Handles unreachable nodes by outputting -1
4. Properly handles weighted directed edges
5. Outputs distances in order of node numbers (1 to n)