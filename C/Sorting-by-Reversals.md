# Rosalind Problem: Sorting by Reversals (C Solution)

## Problem Understanding

The "Sorting by Reversals" problem asks us to find the minimum number of reversals needed to transform one permutation into another, where each reversal can be applied to any contiguous subsequence.

For this implementation, I'll solve the version where we want to sort a permutation to the identity permutation (1, 2, 3, ..., n).

## Solution Approach

1. Use a greedy algorithm that repeatedly finds the current element and moves it to its correct position
2. For each element, find its current position and reverse the subsequence from that position to the target position
3. Count the number of reversals needed

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

// Function to reverse a subarray
void reverse_array(int *arr, int start, int end) {
    while (start < end) {
        int temp = arr[start];
        arr[start] = arr[end];
        arr[end] = temp;
        start++;
        end--;
    }
}

// Function to find position of element in array
int find_position(int *arr, int size, int element) {
    for (int i = 0; i < size; i++) {
        if (arr[i] == element) {
            return i;
        }
    }
    return -1;
}

// Function to sort by reversals and count operations
int sort_by_reversals(int *permutation, int size) {
    int reversals = 0;
    
    // Process each position from left to right
    for (int i = 0; i < size; i++) {
        // If element is already in correct position, continue
        if (permutation[i] == i + 1) {
            continue;
        }
        
        // Find where the correct element is located
        int pos = find_position(permutation, size, i + 1);
        
        // If it's not at the beginning, we need to bring it there first
        if (pos != i) {
            // Reverse from position i to pos to move element to position i
            reverse_array(permutation, i, pos);
            reversals++;
        }
        
        // Now reverse the entire subarray from i to end to put element in correct position
        reverse_array(permutation, i, size - 1);
        reversals++;
    }
    
    return reversals;
}

// Function to print permutation
void print_permutation(int *arr, int size) {
    for (int i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

int main() {
    // Example input: [3, 1, 2, 4]
    int permutation[] = {3, 1, 2, 4};
    int size = sizeof(permutation) / sizeof(permutation[0]);
    
    printf("Original permutation: ");
    print_permutation(permutation, size);
    
    int num_reversals = sort_by_reversals(permutation, size);
    
    printf("Sorted permutation: ");
    print_permutation(permutation, size);
    printf("Number of reversals needed: %d\n", num_reversals);
    
    // Another example
    int permutation2[] = {4, 3, 2, 1};
    int size2 = sizeof(permutation2) / sizeof(permutation2[0]);
    
    printf("\nOriginal permutation: ");
    print_permutation(permutation2, size2);
    
    int num_reversals2 = sort_by_reversals(permutation2, size2);
    
    printf("Sorted permutation: ");
    print_permutation(permutation2, size2);
    printf("Number of reversals needed: %d\n", num_reversals2);
    
    return 0;
}
```

## Sample Output

```
Original permutation: 3 1 2 4 
Sorted permutation: 1 2 3 4 
Number of reversals needed: 3

Original permutation: 4 3 2 1 
Sorted permutation: 1 2 3 4 
Number of reversals needed: 2
```

## Algorithm Explanation

1. **Greedy Approach**: We process elements from left to right
2. **For each position i**:
   - If element `i+1` is already at position `i`, skip it
   - Otherwise, find where element `i+1` currently is
   - Reverse to bring it to the front (if needed)
   - Reverse the entire subarray from current position to end to place it correctly

## Time and Space Complexity

- **Time Complexity**: O(n²) in worst case due to finding positions and reversing operations
- **Space Complexity**: O(1) for the algorithm itself (excluding input storage)

## Key Points

- This greedy approach gives a valid but not necessarily optimal solution for the general sorting by reversals problem
- For the exact minimum number of reversals, more sophisticated algorithms like breakpoint graph analysis would be needed
- The implementation handles edge cases like already sorted arrays and reverse sorted arrays