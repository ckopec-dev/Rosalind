# Rosalind Problem: Longest Increasing Subsequence (LIS) - C Solution

## Problem Understanding
The Longest Increasing Subsequence problem asks us to find the length of the longest subsequence in a given sequence where the elements are in strictly increasing order.

## Approach
We'll use dynamic programming with binary search for an efficient O(n log n) solution:
1. Maintain an array `tails` where `tails[i]` represents the smallest ending element of all increasing subsequences of length `i+1`
2. For each element, use binary search to find the correct position in `tails`
3. Update the `tails` array accordingly

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

// Binary search function to find the position where element should be inserted
int binarySearch(int arr[], int left, int right, int key) {
    while (left < right) {
        int mid = left + (right - left) / 2;
        if (arr[mid] < key) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return left;
}

// Function to find length of longest increasing subsequence
int lisLength(int arr[], int n) {
    if (n == 0) return 0;
    
    // Array to store the smallest tail element for each length
    int *tails = (int*)malloc(n * sizeof(int));
    int length = 0;
    
    for (int i = 0; i < n; i++) {
        // Binary search for the position where arr[i] should be placed
        int pos = binarySearch(tails, 0, length, arr[i]);
        
        // Update the tails array
        tails[pos] = arr[i];
        
        // If we extended the sequence, increment length
        if (pos == length) {
            length++;
        }
    }
    
    free(tails);
    return length;
}

// Function to reconstruct the actual LIS (optional)
void printLIS(int arr[], int n) {
    if (n == 0) return;
    
    int *tails = (int*)malloc(n * sizeof(int));
    int *parent = (int*)malloc(n * sizeof(int));
    int length = 0;
    
    // Initialize parent array
    for (int i = 0; i < n; i++) {
        parent[i] = -1;
    }
    
    for (int i = 0; i < n; i++) {
        int pos = binarySearch(tails, 0, length, arr[i]);
        
        if (pos > 0) {
            parent[i] = pos - 1;  // Store the index of previous element
        }
        
        tails[pos] = arr[i];
        
        if (pos == length) {
            length++;
        }
    }
    
    // Reconstruct the LIS
    int *lis = (int*)malloc(length * sizeof(int));
    int current = length - 1;
    
    for (int i = n - 1; i >= 0; i--) {
        if (current >= 0 && tails[current] == arr[i]) {
            lis[current] = arr[i];
            current--;
        }
    }
    
    printf("Longest Increasing Subsequence: ");
    for (int i = 0; i < length; i++) {
        printf("%d ", lis[i]);
    }
    printf("\n");
    
    free(tails);
    free(parent);
    free(lis);
}

int main() {
    // Example input
    int arr[] = {10, 9, 2, 5, 3, 7, 101, 18};
    int n = sizeof(arr) / sizeof(arr[0]);
    
    printf("Input array: ");
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
    
    int result = lisLength(arr, n);
    printf("Length of Longest Increasing Subsequence: %d\n", result);
    
    printLIS(arr, n);
    
    // Another example
    int arr2[] = {0, 1, 0, 3, 2, 3};
    int n2 = sizeof(arr2) / sizeof(arr2[0]);
    
    printf("\nInput array: ");
    for (int i = 0; i < n2; i++) {
        printf("%d ", arr2[i]);
    }
    printf("\n");
    
    result = lisLength(arr2, n2);
    printf("Length of Longest Increasing Subsequence: %d\n", result);
    
    printLIS(arr2, n2);
    
    return 0;
}
```

## Algorithm Explanation

1. **Key Insight**: We maintain an array `tails` where `tails[i]` stores the smallest ending element of all increasing subsequences of length `i+1`.

2. **Process**:
   - For each element in the input array, we perform binary search on `tails` to find the correct position
   - If the element is larger than all elements in `tails`, we extend the array (increase length)
   - Otherwise, we replace an element in `tails` to maintain the property

3. **Time Complexity**: O(n log n) where n is the length of input array
4. **Space Complexity**: O(n)

## Sample Output
```
Input array: 10 9 2 5 3 7 101 18 
Length of Longest Increasing Subsequence: 4
Longest Increasing Subsequence: 2 3 7 18 

Input array: 0 1 0 3 2 3 
Length of Longest Increasing Subsequence: 4
Longest Increasing Subsequence: 0 1 2 3 
```

## Test Cases

The solution handles various cases:
- Empty array
- Single element array  
- Array with decreasing elements
- Array with increasing elements
- Mixed order arrays

This implementation efficiently solves the Longest Increasing Subsequence problem as required for Rosalind's dataset.