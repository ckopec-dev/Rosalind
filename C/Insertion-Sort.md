# Rosalind Problem: Insertion Sort Solution in C

## Problem Understanding
Insertion sort is a simple sorting algorithm that builds the final sorted array one item at a time. It works by repeatedly taking an element from the unsorted portion and inserting it into its correct position in the sorted portion.

## Solution

```c
#include <stdio.h>
#include <stdlib.h>

// Function to perform insertion sort
void insertion_sort(int arr[], int n) {
    int i, key, j;
    
    // Start from the second element (index 1)
    for (i = 1; i < n; i++) {
        key = arr[i];  // Current element to be inserted
        j = i - 1;     // Index of last element in sorted portion
        
        // Move elements greater than key one position ahead
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        
        // Insert key at its correct position
        arr[j + 1] = key;
    }
}

// Function to print array
void print_array(int arr[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

int main() {
    // Read input from stdin
    int n;
    scanf("%d", &n);
    
    // Allocate memory for array
    int *arr = (int*)malloc(n * sizeof(int));
    
    // Read array elements
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    // Perform insertion sort
    insertion_sort(arr, n);
    
    // Print sorted array
    print_array(arr, n);
    
    // Free allocated memory
    free(arr);
    
    return 0;
}
```

## How it works:

1. **Input Reading**: Read the size of array `n` and then read `n` integers into the array
2. **Insertion Sort Algorithm**:
   - Start from the second element (index 1) since a single element is already "sorted"
   - Take current element as `key`
   - Compare `key` with elements in the sorted portion (to its left)
   - Shift larger elements one position to the right
   - Insert `key` at its correct position
3. **Output**: Print the sorted array

## Time and Space Complexity:
- **Time Complexity**: O(n²) worst case, O(n) best case
- **Space Complexity**: O(1) - sorts in place

## Example Input/Output:
```
Input:
6
6 10 4 5 1 2

Output:
1 2 4 5 6 10
```

This solution handles the standard insertion sort problem and can be compiled and run with any C compiler.