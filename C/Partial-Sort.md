# Rosalind Problem: Partial Sort

## Problem Description
Given: A collection of n (n≤10) integers a_1,a_2,…,a_n and an array of k (k≤10) indices i_1,i_2,…,i_k.

Return: The array A′ such that A′[i] = A[i] for each i≠i_j and A′[i_j] = A[i_j] for each j. In other words, sort the elements at the given indices while keeping all other elements in their original positions.

## Solution

```c
#include <stdio.h>
#include <stdlib.h>

// Comparison function for qsort
int compare(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

int main() {
    int n, k;
    
    // Read the number of elements
    scanf("%d", &n);
    
    // Read the array elements
    int *A = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        scanf("%d", &A[i]);
    }
    
    // Read the number of indices to sort
    scanf("%d", &k);
    
    // Read the indices
    int *indices = (int*)malloc(k * sizeof(int));
    for (int i = 0; i < k; i++) {
        scanf("%d", &indices[i]);
    }
    
    // Create a copy of array A to sort only specified elements
    int *temp = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        temp[i] = A[i];
    }
    
    // Extract elements at specified indices and sort them
    int *sorted_elements = (int*)malloc(k * sizeof(int));
    for (int i = 0; i < k; i++) {
        sorted_elements[i] = temp[indices[i]];
    }
    
    qsort(sorted_elements, k, sizeof(int), compare);
    
    // Put sorted elements back into their positions
    int element_index = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            if (i == indices[j]) {
                A[i] = sorted_elements[element_index];
                element_index++;
                break;
            }
        }
    }
    
    // Print the result
    for (int i = 0; i < n; i++) {
        printf("%d ", A[i]);
    }
    printf("\n");
    
    // Free allocated memory
    free(A);
    free(indices);
    free(temp);
    free(sorted_elements);
    
    return 0;
}
```

## Explanation

1. **Input Reading**: 
   - Read the size of array `n`
   - Read `n` integers into array `A`
   - Read the number of indices `k` to sort
   - Read `k` indices

2. **Sorting Process**:
   - Create a temporary copy of the original array
   - Extract elements at the specified indices
   - Sort these extracted elements using `qsort`
   - Place the sorted elements back into their respective positions in the original array

3. **Output**: Print the modified array where only elements at specified indices are sorted.

## Example Input/Output

**Input:**
```
10
5 2 8 1 9 3 7 4 6 0
4
2 5 7 9
```

**Output:**
```
5 2 0 1 9 3 7 4 6 8 
```

The elements at indices 2, 5, 7, and 9 (which are 8, 3, 4, 0) are sorted to become 0, 3, 4, 8 while other elements remain in their original positions.