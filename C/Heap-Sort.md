# Rosalind Problem: Heap Sort Solution in C

## Problem Understanding

Heap sort is a comparison-based sorting algorithm that uses a binary heap data structure. The algorithm works by:
1. Building a max heap from the input array
2. Repeatedly extracting the maximum element and rebuilding the heap

## Solution Approach

I'll implement heap sort using the standard heapify procedure to maintain the max heap property.

```c
#include <stdio.h>
#include <stdlib.h>

// Function to maintain the max heap property
void heapify(int arr[], int n, int i) {
    int largest = i;        // Initialize largest as root
    int left = 2 * i + 1;   // Left child
    int right = 2 * i + 2;  // Right child
    
    // If left child exists and is greater than root
    if (left < n && arr[left] > arr[largest])
        largest = left;
    
    // If right child exists and is greater than largest so far
    if (right < n && arr[right] > arr[largest])
        largest = right;
    
    // If largest is not root
    if (largest != i) {
        // Swap
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;
        
        // Recursively heapify the affected sub-tree
        heapify(arr, n, largest);
    }
}

// Main heap sort function
void heapSort(int arr[], int n) {
    // Build max heap
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);
    
    // Extract elements from heap one by one
    for (int i = n - 1; i > 0; i--) {
        // Move current root to end
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;
        
        // Call heapify on the reduced heap
        heapify(arr, i, 0);
    }
}

// Function to print array
void printArray(int arr[], int n) {
    for (int i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

int main() {
    // Read input from stdin
    int n;
    scanf("%d", &n);
    
    int *arr = (int*)malloc(n * sizeof(int));
    
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    heapSort(arr, n);
    
    printArray(arr, n);
    
    free(arr);
    return 0;
}
```

## Algorithm Explanation

1. **Build Max Heap**: Starting from the last non-leaf node (index `n/2 - 1`), we heapify each node to create a max heap
2. **Sort Process**: 
   - Extract the maximum element (root) and place it at the end
   - Reduce heap size by 1
   - Heapify the root to restore max heap property
   - Repeat until heap is empty

## Time and Space Complexity

- **Time Complexity**: O(n log n) for both average and worst cases
- **Space Complexity**: O(1) as it's an in-place sorting algorithm

## Sample Input/Output

**Input:**
```
5
4 10 3 5 1
```

**Output:**
```
1 3 4 5 10
```

The solution correctly handles the heap sort algorithm and produces the sorted output as required by Rosalind's problem specifications.