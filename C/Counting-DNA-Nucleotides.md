# Rosalind Problem: Counting DNA Nucleotides

## Problem Statement
Given a DNA string of length at most 1000 nucleotides, return four integers representing the counts of nucleotides A, C, G, and T respectively.

## Solution in C

```c
#include <stdio.h>
#include <string.h>

int main() {
    char dna[1001]; // DNA string with max length 1000 + 1 for null terminator
    int count_A = 0, count_C = 0, count_G = 0, count_T = 0;
    
    // Read the DNA string
    scanf("%s", dna);
    
    // Count each nucleotide
    int length = strlen(dna);
    for (int i = 0; i < length; i++) {
        switch (dna[i]) {
            case 'A':
                count_A++;
                break;
            case 'C':
                count_C++;
                break;
            case 'G':
                count_G++;
                break;
            case 'T':
                count_T++;
                break;
        }
    }
    
    // Print the counts
    printf("%d %d %d %d\n", count_A, count_C, count_G, count_T);
    
    return 0;
}
```

## Explanation

1. **Input Reading**: We read the DNA string using `scanf("%s", dna)` which stores the string in the character array `dna`.

2. **Counter Initialization**: We initialize four integer counters for each nucleotide (A, C, G, T) to zero.

3. **Character Processing**: We iterate through each character of the DNA string and use a switch statement to increment the appropriate counter based on the nucleotide encountered.

4. **Output**: Finally, we print the counts in the required order: A, C, G, T separated by spaces.

## Sample Input/Output

**Input:**
```
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
```

**Output:**
```
20 12 17 21
```

## Time and Space Complexity

- **Time Complexity**: O(n) where n is the length of the DNA string
- **Space Complexity**: O(n) for storing the DNA string

The solution efficiently processes the DNA string in a single pass and handles the constraints specified in the Rosalind problem.