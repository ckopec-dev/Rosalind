# Rosalind Problem: Counting Inversions in C

## Problem Understanding

An inversion in an array is a pair of indices (i, j) such that i < j but arr[i] > arr[j]. We need to count the total number of such inversions in the given array.

## Approach

I'll use the merge sort algorithm to efficiently count inversions in O(n log n) time complexity. During the merge process, when we merge two sorted subarrays, if an element from the right subarray is smaller than an element from the left subarray, it forms inversions with all remaining elements in the left subarray.

## Solution

```c
#include <stdio.h>
#include <stdlib.h>

// Function to merge two sorted arrays and count inversions
long long merge(int *arr, int *temp, int left, int mid, int right) {
    int i, j, k;
    long long invCount = 0;
    
    // Copy data to temp array
    for (i = left; i <= right; i++) {
        temp[i] = arr[i];
    }
    
    i = left;      // Initial index of left subarray
    j = mid + 1;   // Initial index of right subarray
    k = left;      // Initial index of merged array
    
    // Merge the two subarrays while counting inversions
    while (i <= mid && j <= right) {
        if (temp[i] <= temp[j]) {
            arr[k++] = temp[i++];
        } else {
            arr[k++] = temp[j++];
            // All elements from i to mid are greater than temp[j]
            // So we have (mid - i + 1) inversions
            invCount += (mid - i + 1);
        }
    }
    
    // Copy remaining elements
    while (i <= mid) {
        arr[k++] = temp[i++];
    }
    
    while (j <= right) {
        arr[k++] = temp[j++];
    }
    
    return invCount;
}

// Recursive function to count inversions using merge sort
long long mergeSort(int *arr, int *temp, int left, int right) {
    long long invCount = 0;
    
    if (left < right) {
        int mid = left + (right - left) / 2;
        
        invCount += mergeSort(arr, temp, left, mid);
        invCount += mergeSort(arr, temp, mid + 1, right);
        invCount += merge(arr, temp, left, mid, right);
    }
    
    return invCount;
}

// Main function to count inversions
long long countInversions(int *arr, int n) {
    int *temp = (int *)malloc(n * sizeof(int));
    if (temp == NULL) {
        printf("Memory allocation failed\n");
        return -1;
    }
    
    // Make a copy of the original array
    int *arrCopy = (int *)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        arrCopy[i] = arr[i];
    }
    
    long long result = mergeSort(arrCopy, temp, 0, n - 1);
    
    free(temp);
    free(arrCopy);
    
    return result;
}

int main() {
    // Read input
    int n;
    scanf("%d", &n);
    
    int *arr = (int *)malloc(n * sizeof(int));
    if (arr == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }
    
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    // Count inversions
    long long inversions = countInversions(arr, n);
    printf("%lld\n", inversions);
    
    free(arr);
    return 0;
}
```

## Explanation

1. **Merge Sort Approach**: The solution uses a modified merge sort algorithm to count inversions efficiently.

2. **Key Insight**: During the merge step, when we take an element from the right subarray (temp[j]), it means this element is smaller than all remaining elements in the left subarray (temp[i] to temp[mid]). This gives us `(mid - i + 1)` inversions.

3. **Time Complexity**: O(n log n) - same as merge sort
4. **Space Complexity**: O(n) for the temporary array

## Sample Input/Output

**Input:**
```
5
2 3 8 6 1
```

**Output:**
```
5
```

**Explanation**: The inversions are (2,1), (3,1), (8,6), (8,1), (6,1) = 5 inversions.

## How it Works

1. Recursively divide the array into two halves
2. Count inversions in each half
3. Count inversions between the two halves during merge process
4. Return total count of inversions

The algorithm correctly handles all edge cases including already sorted arrays (0 inversions) and reverse sorted arrays (maximum inversions).