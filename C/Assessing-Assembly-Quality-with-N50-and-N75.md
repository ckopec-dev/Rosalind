# Rosalind Problem: Assessing Assembly Quality with N50 and N75

## Problem Understanding

The N50 and N75 statistics are used to assess the quality of genome assemblies:
- **N50**: The smallest contig length such that the sum of contig lengths ≥ 50% of the total assembly
- **N75**: The smallest contig length such that the sum of contig lengths ≥ 75% of the total assembly

## Solution Approach

1. Read contig lengths from input
2. Sort lengths in descending order
3. Calculate cumulative sum
4. Find N50 and N75 values

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

// Comparison function for qsort (descending order)
int compare_desc(const void *a, const void *b) {
    return (*(long long*)b - *(long long*)a);
}

int main() {
    long long n;
    
    // Read number of contigs
    scanf("%lld", &n);
    
    // Allocate memory for contig lengths
    long long *lengths = (long long*)malloc(n * sizeof(long long));
    
    // Read contig lengths
    for (long long i = 0; i < n; i++) {
        scanf("%lld", &lengths[i]);
    }
    
    // Sort in descending order
    qsort(lengths, n, sizeof(long long), compare_desc);
    
    // Calculate total assembly length
    long long total_length = 0;
    for (long long i = 0; i < n; i++) {
        total_length += lengths[i];
    }
    
    // Find N50 and N75
    long long n50 = 0, n75 = 0;
    long long cumulative_sum = 0;
    
    for (long long i = 0; i < n; i++) {
        cumulative_sum += lengths[i];
        
        if (n50 == 0 && cumulative_sum >= total_length / 2) {
            n50 = lengths[i];
        }
        
        if (n75 == 0 && cumulative_sum >= total_length * 3 / 4) {
            n75 = lengths[i];
        }
        
        // If both found, we can break
        if (n50 != 0 && n75 != 0) {
            break;
        }
    }
    
    // Output results
    printf("%lld %lld\n", n50, n75);
    
    // Free allocated memory
    free(lengths);
    
    return 0;
}
```

## How it Works

1. **Input Reading**: Read the number of contigs and their lengths
2. **Sorting**: Sort contig lengths in descending order using `qsort`
3. **Cumulative Sum**: Iterate through sorted lengths, keeping running sum
4. **N50/N75 Calculation**: 
   - N50: First contig length where cumulative sum ≥ 50% of total
   - N75: First contig length where cumulative sum ≥ 75% of total
5. **Output**: Print both values separated by space

## Example Input/Output

**Input:**
```
5
100 200 300 400 500
```

**Output:**
```
400 400
```

## Time Complexity
- Sorting: O(n log n)
- Finding N50/N75: O(n)
- Overall: O(n log n)

## Space Complexity
- O(n) for storing contig lengths

The solution handles the key requirement of finding the minimum contig length such that the sum of all contigs of that length or longer contains at least 50% (N50) or 75% (N75) of the total assembly.