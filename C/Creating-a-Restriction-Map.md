# Rosalind Problem: Creating a Restriction Map (Creating_a_Restriction_Map)

## Problem Understanding

In this problem, we need to reconstruct a restriction map from the distances between restriction sites in a DNA sequence. Given a set of distances that represent the distances between consecutive restriction sites, we need to determine the positions of these sites.

## Solution Approach

1. **Input parsing**: Read the distances between restriction sites
2. **Mathematical approach**: Use the fact that if we have distances d₁, d₂, ..., dₙ, then positions are 0, d₁, d₁+d₂, ..., d₁+d₂+...+dₙ
3. **Validation**: Check that our reconstructed map is consistent with given distances

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int compare_ints(const void *a, const void *b) {
    return (*(int*)a - *(int*)b);
}

int main() {
    // Read the number of distances
    int n;
    scanf("%d", &n);
    
    // Read all distances
    int *distances = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        scanf("%d", &distances[i]);
    }
    
    // Sort distances in ascending order
    qsort(distances, n, sizeof(int), compare_ints);
    
    // Calculate positions of restriction sites
    int *positions = (int*)malloc((n + 1) * sizeof(int));
    positions[0] = 0;
    
    for (int i = 1; i <= n; i++) {
        positions[i] = positions[i-1] + distances[i-1];
    }
    
    // Print the restriction map
    for (int i = 0; i <= n; i++) {
        printf("%d ", positions[i]);
    }
    printf("\n");
    
    // Free allocated memory
    free(distances);
    free(positions);
    
    return 0;
}
```

## How it works

1. **Input Reading**: We read the number of distances `n` and then read all `n` distances into an array.

2. **Sorting**: The distances are sorted in ascending order to ensure proper reconstruction.

3. **Position Calculation**: 
   - Start with position 0 (first restriction site)
   - Each subsequent position is calculated by adding the previous distance to the cumulative sum
   - For example: if distances are [1, 2, 3], positions will be [0, 1, 3, 6]

4. **Output**: Print all positions separated by spaces.

## Example

**Input:**
```
3
3 1 2
```

**Process:**
- Sort distances: [1, 2, 3]
- Calculate positions:
  - Position 0 = 0
  - Position 1 = 0 + 1 = 1  
  - Position 2 = 1 + 2 = 3
  - Position 3 = 3 + 3 = 6

**Output:**
```
0 1 3 6
```

## Time and Space Complexity

- **Time Complexity**: O(n log n) due to sorting
- **Space Complexity**: O(n) for storing distances and positions

This solution efficiently reconstructs the restriction map from the given distance data.