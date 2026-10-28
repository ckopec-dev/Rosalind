```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to calculate factorial
long long factorial(int n) {
    if (n <= 1) return 1;
    return n * factorial(n - 1);
}

// Function to swap two elements
void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Function to print permutation
void printPermutation(int* arr, int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

// Function to generate all permutations using Heap's algorithm
void generatePermutations(int* arr, int size, int* count, int n) {
    if (size == 1) {
        (*count)++;
        printPermutation(arr, n);
        return;
    }
    
    for (int i = 0; i < size; i++) {
        generatePermutations(arr, size - 1, count, n);
        
        if (size % 2 == 1) {
            swap(&arr[0], &arr[size - 1]);
        } else {
            swap(&arr[i], &arr[size - 1]);
        }
    }
}

int main() {
    int n;
    scanf("%d", &n);
    
    // Calculate total number of permutations
    long long total = factorial(n);
    printf("%lld\n", total);
    
    // Create array with numbers 1 to n
    int* arr = (int*)malloc((n + 1) * sizeof(int));
    for (int i = 1; i <= n; i++) {
        arr[i] = i;
    }
    
    int count = 0;
    generatePermutations(arr, n, &count, n);
    
    free(arr);
    return 0;
}
```

This solution implements the Enumerating Gene Orders problem from Rosalind:

1. **Input**: An integer `n` representing the number of genes
2. **Output**: 
   - First line: Total number of permutations (n!)
   - Following lines: All possible permutations of numbers 1 to n

**Key components**:
- `factorial()`: Calculates n! for the total number of permutations
- `swap()`: Helper function to swap array elements
- `printPermutation()`: Prints a single permutation
- `generatePermutations()`: Uses Heap's algorithm to generate all permutations recursively
- `main()`: Reads input, calculates total permutations, and generates/print all permutations

**Time Complexity**: O(n! × n) - since there are n! permutations and each takes O(n) time to print
**Space Complexity**: O(n) - for the array storage and recursion stack

The algorithm uses Heap's algorithm which is efficient for generating permutations without requiring additional memory for storing intermediate results.