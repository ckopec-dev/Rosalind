# Rosalind Problem 2SUM Solution in C

## Problem Understanding
The 2SUM problem asks us to find two distinct indices in an array where the sum of elements at those indices equals zero.

## Solution

```c
#include <stdio.h>
#include <stdlib.h>

// Function to compare integers for qsort
int compare(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

// Function to find two sum indices
int* two_sum(int arr[], int n) {
    // Create array of pairs (value, original_index)
    typedef struct {
        int value;
        int index;
    } pair;
    
    pair *pairs = malloc(n * sizeof(pair));
    
    // Fill the pairs array
    for (int i = 0; i < n; i++) {
        pairs[i].value = arr[i];
        pairs[i].index = i + 1; // 1-indexed
    }
    
    // Sort by value
    qsort(pairs, n, sizeof(pair), compare);
    
    // Two pointers approach
    int left = 0;
    int right = n - 1;
    int *result = malloc(2 * sizeof(int));
    
    while (left < right) {
        int sum = pairs[left].value + pairs[right].value;
        
        if (sum == 0) {
            // Found the pair
            result[0] = pairs[left].index;
            result[1] = pairs[right].index;
            
            // Ensure result is in ascending order
            if (result[0] > result[1]) {
                int temp = result[0];
                result[0] = result[1];
                result[1] = temp;
            }
            
            free(pairs);
            return result;
        } else if (sum < 0) {
            left++;
        } else {
            right--;
        }
    }
    
    // No solution found
    free(pairs);
    return NULL;
}

int main() {
    int n;
    scanf("%d", &n);
    
    int *arr = malloc(n * sizeof(int));
    
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    int *result = two_sum(arr, n);
    
    if (result != NULL) {
        printf("%d %d\n", result[0], result[1]);
        free(result);
    } else {
        printf("No solution found\n");
    }
    
    free(arr);
    return 0;
}
```

## Explanation

1. **Data Structure**: We create a structure to hold both the value and original index of each element.

2. **Sorting**: We sort the array while keeping track of original indices.

3. **Two Pointers Technique**: 
   - Use left pointer at start and right pointer at end
   - If sum equals 0, we found our answer
   - If sum is less than 0, move left pointer right (increase sum)
   - If sum is greater than 0, move right pointer left (decrease sum)

4. **Output**: Return indices in ascending order as required.

## Time and Space Complexity

- **Time Complexity**: O(n log n) due to sorting
- **Space Complexity**: O(n) for the additional pairs array

## Sample Input/Output

**Input:**
```
6
2 -3 4 -1 5 -2
```

**Output:**
```
2 6
```

This represents indices (1-indexed) where elements -3 and -2 sum to zero.