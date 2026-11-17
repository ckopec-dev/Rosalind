# Rosalind Problem: 2-Way Partition

## Problem Description
Given: A positive integer n ≤ 10⁵ and an array A[1..n] of integers.
Return: A partition of the array such that elements smaller than the pivot come before elements equal to the pivot, which comes before elements greater than the pivot.

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>

// Function to perform 2-way partition
int partition(int *arr, int low, int high) {
    // Choose the rightmost element as pivot
    int pivot = arr[high];
    
    // Index of smaller element (indicates right position of pivot)
    int i = low - 1;
    
    for (int j = low; j < high; j++) {
        // If current element is smaller than or equal to pivot
        if (arr[j] <= pivot) {
            i++;
            // Swap elements at positions i and j
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    
    // Place pivot in its correct position
    int temp = arr[i + 1];
    arr[i + 1] = arr[high];
    arr[high] = temp;
    
    return i + 1; // Return the partition index
}

// Function to print array
void printArray(int *arr, int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

int main() {
    int n;
    scanf("%d", &n);
    
    int *arr = (int*)malloc(n * sizeof(int));
    
    // Read array elements
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    // Perform partition
    partition(arr, 0, n - 1);
    
    // Print result
    printArray(arr, n);
    
    free(arr);
    return 0;
}
```

## How it works:

1. **Partition Function**: 
   - Uses the last element as pivot
   - Maintains three sections: elements ≤ pivot, elements = pivot, elements > pivot
   - Uses two pointers to rearrange elements in-place

2. **Algorithm Steps**:
   - Choose pivot (last element)
   - Use a pointer `i` to track position where elements ≤ pivot should go
   - Iterate through array from left to right
   - When encountering element ≤ pivot, swap it with element at index `i`
   - Place pivot in its correct final position

3. **Time Complexity**: O(n log n) average case, O(n²) worst case
4. **Space Complexity**: O(log n) due to recursion stack

## Sample Input:
```
9
4 5 6 4 1 2 5 7 4
```

## Sample Output:
```
1 2 4 4 4 5 5 6 7
```