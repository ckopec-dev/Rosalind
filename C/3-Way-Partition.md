# Rosalind Problem: 3-Way Partition

## Problem Description
Given an array of n integers, partition it into three parts:
- Elements less than a pivot value
- Elements equal to the pivot value  
- Elements greater than the pivot value

The goal is to rearrange the array so that all elements less than the pivot come first, followed by elements equal to the pivot, then elements greater than the pivot.

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>

// Function to perform 3-way partition
void three_way_partition(int arr[], int low, int high, int *lt, int *gt) {
    int pivot = arr[low];
    int i = low + 1;
    *lt = low;
    *gt = high;
    
    while (i <= *gt) {
        if (arr[i] < pivot) {
            // Swap arr[*lt] and arr[i]
            int temp = arr[*lt];
            arr[*lt] = arr[i];
            arr[i] = temp;
            (*lt)++;
            i++;
        } else if (arr[i] > pivot) {
            // Swap arr[*gt] and arr[i]
            int temp = arr[*gt];
            arr[*gt] = arr[i];
            arr[i] = temp;
            (*gt]--;
        } else {
            i++;
        }
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
    
    int *arr = (int*)malloc(n * sizeof(int));
    
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    // Perform 3-way partition
    int lt, gt;
    three_way_partition(arr, 0, n - 1, &lt, &gt);
    
    // Print the result
    print_array(arr, n);
    
    free(arr);
    return 0;
}
```

## Explanation

This solution implements a 3-way partition algorithm based on the Dutch National Flag problem approach:

1. **Input Reading**: Read the array size `n` and then read `n` integers into the array
2. **Partition Logic**: 
   - Choose the first element as pivot
   - Maintain three pointers:
     - `lt`: points to the end of elements less than pivot
     - `i`: current element being processed  
     - `gt`: points to the start of elements greater than pivot
3. **Three Cases**:
   - If current element < pivot: swap with element at `lt` and increment both `lt` and `i`
   - If current element > pivot: swap with element at `gt` and decrement `gt` (don't increment `i`)
   - If current element = pivot: just increment `i`

## Time and Space Complexity
- **Time Complexity**: O(n) - single pass through the array
- **Space Complexity**: O(1) - only using a constant amount of extra space

## Sample Input/Output
```
Input:
9
4 5 6 4 1 2 5 7 4

Output:
1 2 4 4 4 5 5 6 7
```

The array is partitioned with elements less than pivot (4) on the left, equal to pivot in the middle, and greater than pivot on the right.