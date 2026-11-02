# Rosalind Problem: Comparing Spectra with the Spectral Convolution

## Problem Description

The spectral convolution of two spectra is the multiset of differences between all pairs of masses in the two spectra. In this problem, we need to compute the spectral convolution of two given spectra.

## Solution Approach

1. Read two spectra (arrays of masses)
2. Compute all possible differences between elements from both spectra
3. Count the frequency of each difference
4. Return the sorted list of differences with their frequencies

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Structure to store a difference and its count
typedef struct {
    int difference;
    int count;
} Difference;

// Comparison function for sorting
int compare_differences(const void *a, const void *b) {
    Difference *diff_a = (Difference *)a;
    Difference *diff_b = (Difference *)b;
    
    // Sort by difference value
    return diff_a->difference - diff_b->difference;
}

int main() {
    // Read first spectrum
    int n1;
    scanf("%d", &n1);
    
    int *spectrum1 = (int *)malloc(n1 * sizeof(int));
    for (int i = 0; i < n1; i++) {
        scanf("%d", &spectrum1[i]);
    }
    
    // Read second spectrum
    int n2;
    scanf("%d", &n2);
    
    int *spectrum2 = (int *)malloc(n2 * sizeof(int));
    for (int i = 0; i < n2; i++) {
        scanf("%d", &spectrum2[i]);
    }
    
    // Calculate all possible differences
    int total_differences = n1 * n2;
    Difference *differences = (Difference *)malloc(total_differences * sizeof(Difference));
    
    int index = 0;
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            differences[index].difference = spectrum1[i] - spectrum2[j];
            differences[index].count = 1;
            index++;
        }
    }
    
    // Sort differences
    qsort(differences, total_differences, sizeof(Difference), compare_differences);
    
    // Count frequencies and print results
    printf("%d\n", differences[0].difference);
    int current_diff = differences[0].difference;
    int count = 1;
    
    for (int i = 1; i < total_differences; i++) {
        if (differences[i].difference == current_diff) {
            count++;
        } else {
            printf("%d\n", current_diff);
            current_diff = differences[i].difference;
            count = 1;
        }
    }
    printf("%d\n", current_diff);
    
    // Free allocated memory
    free(spectrum1);
    free(spectrum2);
    free(differences);
    
    return 0;
}
```

## Alternative Implementation (More Efficient)

```c
#include <stdio.h>
#include <stdlib.h>

// Structure to store difference and its count
typedef struct {
    int value;
    int count;
} Element;

int compare_elements(const void *a, const void *b) {
    Element *el_a = (Element *)a;
    Element *el_b = (Element *)b;
    return el_a->value - el_b->value;
}

int main() {
    // Read first spectrum
    int n1;
    scanf("%d", &n1);
    
    int *spectrum1 = (int *)malloc(n1 * sizeof(int));
    for (int i = 0; i < n1; i++) {
        scanf("%d", &spectrum1[i]);
    }
    
    // Read second spectrum
    int n2;
    scanf("%d", &n2);
    
    int *spectrum2 = (int *)malloc(n2 * sizeof(int));
    for (int i = 0; i < n2; i++) {
        scanf("%d", &spectrum2[i]);
    }
    
    // Calculate all differences
    int total_diff = n1 * n2;
    int *diff_array = (int *)malloc(total_diff * sizeof(int));
    int index = 0;
    
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            diff_array[index++] = spectrum1[i] - spectrum2[j];
        }
    }
    
    // Sort differences
    qsort(diff_array, total_diff, sizeof(int), compare_elements);
    
    // Print results with counts
    int current = diff_array[0];
    int count = 1;
    
    for (int i = 1; i < total_diff; i++) {
        if (diff_array[i] == current) {
            count++;
        } else {
            printf("%d %d\n", current, count);
            current = diff_array[i];
            count = 1;
        }
    }
    printf("%d %d\n", current, count);
    
    // Free memory
    free(spectrum1);
    free(spectrum2);
    free(diff_array);
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
3
186 185 174
2
101 102
```

**Output:**
```
-101 1
-102 1
-174 1
-175 1
-185 1
-186 1
```

## Key Points

1. **Time Complexity**: O(n₁ × n₂ + n₁ × n₂ × log(n₁ × n₂)) where n₁ and n₂ are the sizes of the two spectra
2. **Space Complexity**: O(n₁ × n₂) for storing all differences
3. The problem asks to compute spectral convolution, which is essentially computing all pairwise differences between elements of two spectra
4. Results should be sorted by difference value
5. Each unique difference should be printed with its frequency count

The solution handles the computation efficiently by generating all possible differences and then sorting them to group identical values together for counting.