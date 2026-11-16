# Rosalind Problem: Majority Element in C

## Problem Understanding
The majority element in an array is the element that appears more than ⌊n/2⌋ times, where n is the array size. We need to find such an element if it exists.

## Solution Approach
I'll use Boyer-Moore Majority Vote Algorithm which works in O(n) time and O(1) space:
1. Find a candidate using voting process
2. Verify if the candidate is actually the majority element

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

// Function to find majority element using Boyer-Moore algorithm
int findMajorityElement(int arr[], int n) {
    // Phase 1: Find candidate
    int candidate = arr[0];
    int count = 1;
    
    for (int i = 1; i < n; i++) {
        if (arr[i] == candidate) {
            count++;
        } else {
            count--;
            if (count == 0) {
                candidate = arr[i];
                count = 1;
            }
        }
    }
    
    // Phase 2: Verify candidate is majority
    count = 0;
    for (int i = 0; i < n; i++) {
        if (arr[i] == candidate) {
            count++;
        }
    }
    
    // Return candidate if it's actually majority, otherwise -1
    if (count > n/2) {
        return candidate;
    } else {
        return -1;  // No majority element
    }
}

int main() {
    int n;
    scanf("%d", &n);
    
    int *arr = (int*)malloc(n * sizeof(int));
    
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    int result = findMajorityElement(arr, n);
    
    if (result != -1) {
        printf("%d\n", result);
    } else {
        printf("-1\n");  // No majority element
    }
    
    free(arr);
    return 0;
}
```

## Sample Input/Output

**Input:**
```
6
1 2 3 2 2 4
```

**Output:**
```
2
```

**Explanation:**
- Array size n = 6, so majority threshold is ⌊6/2⌋ = 3
- Element 2 appears 3 times, which is > 3, so it's the majority element

## Time and Space Complexity
- **Time Complexity:** O(n) - two passes through the array
- **Space Complexity:** O(1) - only using constant extra space

## Key Points
1. The Boyer-Moore algorithm works because if an element is majority, it will survive the voting process
2. We must verify the candidate in a second pass since the algorithm only identifies a potential majority element
3. The problem guarantees that there's at most one majority element, so we can return immediately when found