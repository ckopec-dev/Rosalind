# Rosalind Problem: Binary Search in C

## Problem Description
Given: A positive integer n ≤ 10^5 and a sorted array A[1..n] of integers from -10^5 to 10^5, followed by a positive integer k ≤ 10^5 and k integers from -10^5 to 10^5.

Return: For each integer k in the second list, output the index (1-based) of the integer in the array A. If the integer is not found, output -1.

## Solution

```c
#include <stdio.h>

// Binary search function that returns 1-based index or -1 if not found
int binary_search(int arr[], int left, int right, int target) {
    while (left <= right) {
        int mid = left + (right - left) / 2;
        
        if (arr[mid] == target) {
            return mid + 1; // Return 1-based index
        }
        else if (arr[mid] < target) {
            left = mid + 1;
        }
        else {
            right = mid - 1;
        }
    }
    return -1; // Not found
}

int main() {
    int n, k;
    
    // Read the size of array
    scanf("%d", &n);
    
    // Read the sorted array
    int arr[n];
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    // Read the number of elements to search
    scanf("%d", &k);
    
    // Read and search each element
    for (int i = 0; i < k; i++) {
        int target;
        scanf("%d", &target);
        
        int result = binary_search(arr, 0, n - 1, target);
        printf("%d ", result);
    }
    
    printf("\n");
    return 0;
}
```

## Explanation

### Key Points:
1. **Binary Search Implementation**: The `binary_search` function implements the classic binary search algorithm with O(log n) time complexity.

2. **Index Handling**: 
   - Array indices are 0-based in C
   - Problem requires 1-based indexing, so we add 1 to the returned index
   - If element not found, return -1

3. **Algorithm Steps**:
   - Calculate middle index using `mid = left + (right - left) / 2` to avoid overflow
   - Compare target with middle element
   - Adjust search range based on comparison
   - Return index when found, or -1 if not found

### Input/Output Example:
```
Input:
5
1 2 3 4 5
2
3 6

Output:
3 -1
```

### Time Complexity: O(k log n)
- Where k is the number of elements to search
- Each binary search takes O(log n) time

### Space Complexity: O(n)
- For storing the array of size n