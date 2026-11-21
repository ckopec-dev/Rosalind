# Rosalind Problem: Median (MED)

## Problem Description
Given a collection of n (n≤10^5) positive integers, return the median of the collection.

The median is the middle value when the collection is sorted. If there are an even number of elements, the median is the average of the two middle elements.

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>

// Comparison function for qsort
int compare(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

int main() {
    int n;
    
    // Read the number of integers
    scanf("%d", &n);
    
    // Allocate memory for array
    int *arr = (int*)malloc(n * sizeof(int));
    
    // Read all integers
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    // Sort the array
    qsort(arr, n, sizeof(int), compare);
    
    // Calculate median
    double median;
    if (n % 2 == 1) {
        // Odd number of elements - median is the middle element
        median = arr[n/2];
    } else {
        // Even number of elements - median is average of two middle elements
        median = (arr[n/2 - 1] + arr[n/2]) / 2.0;
    }
    
    // Print the median
    printf("%.1f\n", median);
    
    // Free allocated memory
    free(arr);
    
    return 0;
}
```

## Explanation

1. **Input Reading**: First read the number of integers `n`, then read all `n` integers into an array.

2. **Sorting**: Use `qsort()` function with a custom comparison function to sort the array in ascending order.

3. **Median Calculation**:
   - If `n` is odd, the median is the element at index `n/2`
   - If `n` is even, the median is the average of elements at indices `n/2-1` and `n/2`

4. **Output**: Print the median with one decimal place as required by Rosalind format.

## Sample Input/Output

**Input:**
```
5
3 1 4 1 5
```

**Output:**
```
3.0
```

**Explanation**: After sorting [1, 1, 3, 4, 5], the median is the middle element (index 2) which is 3.

## Time and Space Complexity

- **Time Complexity**: O(n log n) due to sorting
- **Space Complexity**: O(n) for storing the array

The solution handles the constraints efficiently and produces the correct median value as required by the Rosalind problem.