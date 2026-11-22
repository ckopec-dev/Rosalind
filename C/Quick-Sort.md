# Rosalind Problem: Quick_Sort Solution in C

## Problem Understanding
The quick sort problem requires implementing the quicksort algorithm to sort an array of integers in ascending order.

## Solution Approach
I'll implement the classic quicksort algorithm using:
- A recursive partitioning approach
- Choosing the first element as pivot
- Partitioning elements smaller and larger than pivot
- Recursively sorting subarrays

## C Implementation

```c
#include <stdio.h>

// Function to swap two elements
void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Partition function - arranges elements around pivot
int partition(int arr[], int low, int high) {
    // Choose first element as pivot
    int pivot = arr[low];
    int i = low + 1;  // Index of smaller element
    int j = high;     // Index of larger element
    
    while (i <= j) {
        // Find element greater than or equal to pivot from left
        while (i <= j && arr[i] <= pivot) {
            i++;
        }
        
        // Find element smaller than or equal to pivot from right
        while (i <= j && arr[j] > pivot) {
            j--;
        }
        
        // Swap if elements found on both sides
        if (i < j) {
            swap(&arr[i], &arr[j]);
            i++;
            j--;
        }
    }
    
    // Place pivot in correct position
    swap(&arr[low], &arr[j]);
    return j;
}

// Quick sort function
void quickSort(int arr[], int low, int high) {
    if (low < high) {
        // Partition the array and get pivot index
        int pi = partition(arr, low, high);
        
        // Recursively sort elements before and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

// Main function to read input and solve the problem
int main() {
    int n;
    
    // Read array size
    scanf("%d", &n);
    
    // Read array elements
    int arr[n];
    for (int i = 0; i < n; i++) {
        scanf("%d", &arr[i]);
    }
    
    // Sort the array using quicksort
    quickSort(arr, 0, n - 1);
    
    // Print sorted array
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
    
    return 0;
}
```

## How it works:

1. **Partition Function**: 
   - Selects the first element as pivot
   - Rearranges array so elements smaller than pivot come before it, and larger elements come after
   - Returns the final position of the pivot

2. **QuickSort Function**:
   - Recursively applies quicksort to subarrays on both sides of the pivot
   - Base case: when low >= high, array is already sorted

3. **Main Function**:
   - Reads input array size and elements
   - Calls quickSort to sort the array
   - Prints the sorted result

## Time Complexity:
- Average case: O(n log n)
- Worst case: O(n²) 
- Best case: O(n log n)

## Space Complexity:
- O(log n) due to recursion stack

## Sample Input/Output:

**Input:**
```
10
8 2 4 7 1 3 9 6 5 0
```

**Output:**
```
0 1 2 3 4 5 6 7 8 9 
```