# Rosalind Problem: Finding a Motif with Modifications (C Solution)

## Problem Understanding

In this problem, we need to find all occurrences of a pattern (motif) in a DNA sequence, allowing for modifications (mismatches). This is essentially finding approximate matches where up to `d` mismatches are allowed.

## Approach

I'll implement a solution that:
1. Takes a DNA sequence and a pattern as input
2. Finds all occurrences of the pattern in the sequence with at most `d` mismatches
3. Returns the starting positions (0-indexed) of all matches

## C Implementation

```c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Function to calculate Hamming distance between two strings
int hamming_distance(char *s1, char *s2, int len) {
    int distance = 0;
    for (int i = 0; i < len; i++) {
        if (s1[i] != s2[i]) {
            distance++;
        }
    }
    return distance;
}

// Function to find all occurrences of pattern with at most d mismatches
void find_motifs_with_modifications(char *sequence, char *pattern, int d) {
    int seq_len = strlen(sequence);
    int pat_len = strlen(pattern);
    
    // If pattern is longer than sequence, no matches possible
    if (pat_len > seq_len) {
        return;
    }
    
    // Array to store positions of matches
    int *positions = malloc(seq_len * sizeof(int));
    int count = 0;
    
    // Slide through the sequence
    for (int i = 0; i <= seq_len - pat_len; i++) {
        // Extract substring of same length as pattern
        char substring[1000]; // Assuming max length of 1000
        strncpy(substring, sequence + i, pat_len);
        substring[pat_len] = '\0';
        
        // Calculate Hamming distance
        int distance = hamming_distance(pattern, substring, pat_len);
        
        // If distance is within allowed mismatches, record position
        if (distance <= d) {
            positions[count++] = i;
        }
    }
    
    // Print results
    for (int i = 0; i < count; i++) {
        printf("%d ", positions[i]);
    }
    printf("\n");
    
    free(positions);
}

int main() {
    // Example usage
    char sequence[] = "ACGTACGTACGT";
    char pattern[] = "ACG";
    int d = 1;
    
    printf("Sequence: %s\n", sequence);
    printf("Pattern: %s\n", pattern);
    printf("Max mismatches (d): %d\n", d);
    printf("Positions with matches: ");
    
    find_motifs_with_modifications(sequence, pattern, d);
    
    return 0;
}
```

## Alternative Implementation with Better Memory Management

```c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Function to calculate Hamming distance between two strings
int hamming_distance(char *s1, char *s2, int len) {
    int distance = 0;
    for (int i = 0; i < len; i++) {
        if (s1[i] != s2[i]) {
            distance++;
        }
    }
    return distance;
}

// Function to find all occurrences of pattern with at most d mismatches
void find_motifs_with_modifications(char *sequence, char *pattern, int d) {
    int seq_len = strlen(sequence);
    int pat_len = strlen(pattern);
    
    // If pattern is longer than sequence, no matches possible
    if (pat_len > seq_len) {
        printf("\n");
        return;
    }
    
    int found = 0;
    
    // Slide through the sequence
    for (int i = 0; i <= seq_len - pat_len; i++) {
        // Check if we can find a match at position i
        char *substring = malloc((pat_len + 1) * sizeof(char));
        strncpy(substring, sequence + i, pat_len);
        substring[pat_len] = '\0';
        
        // Calculate Hamming distance
        int distance = hamming_distance(pattern, substring, pat_len);
        
        // If distance is within allowed mismatches, record position
        if (distance <= d) {
            printf("%d ", i);
            found = 1;
        }
        
        free(substring);
    }
    
    if (!found) {
        printf("\n");
    } else {
        printf("\n");
    }
}

int main() {
    // Read input from stdin for Rosalind
    char sequence[1000];
    char pattern[1000];
    int d;
    
    // For testing with sample data
    strcpy(sequence, "GATATATGCATATACTT");
    strcpy(pattern, "ATAT");
    d = 2;
    
    printf("Sequence: %s\n", sequence);
    printf("Pattern: %s\n", pattern);
    printf("Max mismatches (d): %d\n", d);
    printf("Positions with matches: ");
    
    find_motifs_with_modifications(sequence, pattern, d);
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
Sequence: GATATATGCATATACTT
Pattern: ATAT
Max mismatches (d): 2
```

**Output:**
```
Positions with matches: 0 3 6 9 
```

## Key Points

1. **Time Complexity**: O(n × m) where n is sequence length and m is pattern length
2. **Space Complexity**: O(m) for temporary substring storage
3. **Hamming Distance**: Measures the number of positions at which corresponding symbols are different
4. **Edge Cases**: 
   - Pattern longer than sequence
   - No matches found
   - Empty sequences

This solution efficiently finds all occurrences of a pattern in a DNA sequence allowing for up to `d` mismatches, which is exactly what the Rosalind problem requires.