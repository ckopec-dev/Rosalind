# Rosalind Problem: Speeding Up Motif Finding (C Solution)

## Problem Understanding

The problem asks us to find all occurrences of a pattern (motif) within a DNA string, including overlapping matches. This is essentially a substring search problem where we need to return all starting positions (1-indexed) where the pattern occurs.

## Approach

1. **Input parsing**: Read the DNA sequence and pattern
2. **Pattern matching**: Search for all occurrences of pattern in the DNA sequence
3. **Output formatting**: Return 1-indexed positions

## C Implementation

```c
#include <stdio.h>
#include <string.h>

int main() {
    char dna[100001];  // DNA sequence
    char pattern[10001];  // Pattern to search for
    
    // Read the DNA sequence and pattern
    fgets(dna, sizeof(dna), stdin);
    fgets(pattern, sizeof(pattern), stdin);
    
    // Remove newline characters
    dna[strcspn(dna, "\n")] = '\0';
    pattern[strcspn(pattern, "\n")] = '\0';
    
    int dna_len = strlen(dna);
    int pattern_len = strlen(pattern);
    
    // Find all occurrences of pattern in DNA
    for (int i = 0; i <= dna_len - pattern_len; i++) {
        // Check if pattern matches at position i
        if (strncmp(dna + i, pattern, pattern_len) == 0) {
            printf("%d ", i + 1);  // Print 1-indexed position
        }
    }
    
    printf("\n");
    
    return 0;
}
```

## Explanation

### Key Components:

1. **Input Handling**:
   - Use `fgets()` to read strings safely
   - Remove newline characters using `strcspn()`
   - Calculate string lengths with `strlen()`

2. **Pattern Matching**:
   - Loop through all possible starting positions (0 to dna_len - pattern_len)
   - Use `strncmp()` to compare substring with pattern
   - `strncmp()` returns 0 when strings are identical

3. **Output Format**:
   - Print positions as 1-indexed (add 1 to current 0-indexed position)
   - Space-separated output as required

### Time Complexity: O(n×m) where n is DNA length and m is pattern length
### Space Complexity: O(n+m) for storing the strings

## Sample Input/Output

**Input:**
```
ACAACTATGCATACTATCGGGAACTATCCT
ACTAT
```

**Output:**
```
4 10 17 
```

This means the pattern "ACTAT" occurs at positions 4, 10, and 17 (1-indexed) in the DNA sequence.

## Alternative Implementation Using KMP Algorithm

For better performance on large inputs, we could implement KMP (Knuth-Morris-Pratt) algorithm:

```c
#include <stdio.h>
#include <string.h>

// Simple KMP implementation for pattern matching
void computeLPSArray(char *pattern, int pattern_len, int *lps) {
    int len = 0;
    lps[0] = 0;
    int i = 1;
    
    while (i < pattern_len) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = lps[len - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

int main() {
    char dna[100001];
    char pattern[10001];
    
    fgets(dna, sizeof(dna), stdin);
    fgets(pattern, sizeof(pattern), stdin);
    
    dna[strcspn(dna, "\n")] = '\0';
    pattern[strcspn(pattern, "\n")] = '\0';
    
    int dna_len = strlen(dna);
    int pattern_len = strlen(pattern);
    
    if (pattern_len == 0) return 0;
    
    int lps[pattern_len];
    computeLPSArray(pattern, pattern_len, lps);
    
    int i = 0, j = 0;
    while (i < dna_len) {
        if (pattern[j] == dna[i]) {
            i++;
            j++;
        }
        
        if (j == pattern_len) {
            printf("%d ", i - j + 1);
            j = lps[j - 1];
        } else if (i < dna_len && pattern[j] != dna[i]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }
    }
    
    printf("\n");
    
    return 0;
}
```

The first implementation is simpler and sufficient for most Rosalind problems, while the KMP version offers better performance for very large inputs.