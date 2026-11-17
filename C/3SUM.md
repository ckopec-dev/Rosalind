# Rosalind Problem 3SUM Solution in C

## Problem Understanding
The 3SUM problem asks us to find all unique triplets in an array that sum to zero. Given an array of integers, we need to return all combinations of three numbers that add up to zero.

## Solution Approach
1. Sort the array first
2. For each element, use two pointers technique to find pairs that sum to the negative of that element
3. Skip duplicates to ensure unique triplets

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

// Comparison function for qsort
int compare(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

// Function to find all 3SUM combinations
int** threeSum(int* nums, int numsSize, int* returnSize, int** returnColumnSizes) {
    // Sort the array first
    qsort(nums, numsSize, sizeof(int), compare);
    
    // Allocate memory for result
    int** result = (int**)malloc(sizeof(int*) * numsSize * numsSize);
    *returnColumnSizes = (int*)malloc(sizeof(int) * numsSize * numsSize);
    *returnSize = 0;
    
    // Iterate through each element as the first element of triplet
    for (int i = 0; i < numsSize - 2; i++) {
        // Skip duplicates for the first element
        if (i > 0 && nums[i] == nums[i-1]) {
            continue;
        }
        
        // Two pointers approach
        int left = i + 1;
        int right = numsSize - 1;
        
        while (left < right) {
            int sum = nums[i] + nums[left] + nums[right];
            
            if (sum == 0) {
                // Found a triplet
                result[*returnSize] = (int*)malloc(sizeof(int) * 3);
                result[*returnSize][0] = nums[i];
                result[*returnSize][1] = nums[left];
                result[*returnSize][2] = nums[right];
                (*returnColumnSizes)[*returnSize] = 3;
                (*returnSize)++;
                
                // Skip duplicates for second element
                int leftVal = nums[left];
                while (left < right && nums[left] == leftVal) {
                    left++;
                }
                
                // Skip duplicates for third element
                int rightVal = nums[right];
                while (left < right && nums[right] == rightVal) {
                    right--;
                }
            } else if (sum < 0) {
                // Sum is too small, move left pointer right
                left++;
            } else {
                // Sum is too large, move right pointer left
                right--;
            }
        }
    }
    
    return result;
}

// Helper function to print results
void printResult(int** result, int returnSize, int* returnColumnSizes) {
    for (int i = 0; i < returnSize; i++) {
        printf("[%d,%d,%d]\n", result[i][0], result[i][1], result[i][2]);
    }
}

// Example usage
int main() {
    // Test case from Rosalind
    int nums[] = {-1, 0, 1, 2, -1, -4};
    int numsSize = 6;
    
    int returnSize;
    int* returnColumnSizes;
    
    int** result = threeSum(nums, numsSize, &returnSize, &returnColumnSizes);
    
    printf("3SUM results:\n");
    printResult(result, returnSize, returnColumnSizes);
    
    // Free allocated memory
    for (int i = 0; i < returnSize; i++) {
        free(result[i]);
    }
    free(result);
    free(returnColumnSizes);
    
    return 0;
}
```

## Key Features of the Solution

1. **Sorting**: First sorts the array to enable efficient two-pointer technique
2. **Two-pointer Technique**: For each element, uses left and right pointers to find pairs
3. **Duplicate Handling**: Skips duplicate elements to ensure unique triplets
4. **Memory Management**: Properly allocates and frees memory for results

## Time and Space Complexity

- **Time Complexity**: O(n²) where n is the size of the input array
  - Sorting: O(n log n)
  - Two-pointer search: O(n²)
- **Space Complexity**: O(1) excluding the output space

## Sample Output
For input `[-1, 0, 1, 2, -1, -4]`, the output would be:
```
[-1,-1,2]
[-1,0,1]
```

This solution efficiently handles the 3SUM problem and can be adapted for various test cases including edge cases with many duplicates.