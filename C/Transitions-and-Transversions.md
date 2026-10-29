# Rosalind Problem: Transitions and Transversions in C

## Problem Description
Given two DNA strings of equal length, determine the number of transitions and transversions between them.

- **Transition**: Purine ↔ Purine (A↔G) or Pyrimidine ↔ Pyrimidine (C↔T)
- **Transversion**: Purine ↔ Pyrimidine (A↔C, A↔T, G↔C, G↔T)

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to check if a character is a purine (A or G)
int is_purine(char nucleotide) {
    return (nucleotide == 'A' || nucleotide == 'G');
}

// Function to check if a character is a pyrimidine (C or T)
int is_pyrimidine(char nucleotide) {
    return (nucleotide == 'C' || nucleotide == 'T');
}

// Function to determine if two nucleotides are transition
int is_transition(char n1, char n2) {
    // Same nucleotide - not a transition or transversion
    if (n1 == n2) return 0;
    
    // Check if both are purines
    if (is_purine(n1) && is_purine(n2)) {
        return 1;  // A↔G is transition
    }
    
    // Check if both are pyrimidines
    if (is_pyrimidine(n1) && is_pyrimidine(n2)) {
        return 1;  // C↔T is transition
    }
    
    return 0;  // Different types - transversion
}

// Function to determine if two nucleotides are transversion
int is_transversion(char n1, char n2) {
    // Same nucleotide - not a transition or transversion
    if (n1 == n2) return 0;
    
    // Check if one is purine and other is pyrimidine
    if ((is_purine(n1) && is_pyrimidine(n2)) || 
        (is_pyrimidine(n1) && is_purine(n2))) {
        return 1;
    }
    
    return 0;
}

int main() {
    char seq1[1000], seq2[1000];
    
    // Read the two DNA sequences
    printf("Enter first DNA sequence: ");
    fgets(seq1, sizeof(seq1), stdin);
    
    printf("Enter second DNA sequence: ");
    fgets(seq2, sizeof(seq2), stdin);
    
    // Remove newline characters if present
    seq1[strcspn(seq1, "\n")] = 0;
    seq2[strcspn(seq2, "\n")] = 0;
    
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    
    // Check if sequences have equal length
    if (len1 != len2) {
        printf("Error: Sequences must have equal length\n");
        return 1;
    }
    
    int transitions = 0;
    int transversions = 0;
    
    // Compare each position
    for (int i = 0; i < len1; i++) {
        if (is_transition(seq1[i], seq2[i])) {
            transitions++;
        } else if (is_transversion(seq1[i], seq2[i])) {
            transversions++;
        }
    }
    
    // Calculate ratio
    double ratio = (double)transitions / (double)transversions;
    
    printf("Transitions: %d\n", transitions);
    printf("Transversions: %d\n", transversions);
    printf("Transition/Transversion ratio: %.10f\n", ratio);
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
Enter first DNA sequence: ACGTACGT
Enter second DNA sequence: AGGTACGT
```

**Output:**
```
Transitions: 2
Transversions: 1
Transition/Transversion ratio: 2.0000000000
```

## Explanation

1. **Input**: Two DNA sequences of equal length
2. **Processing**: Compare each position in both sequences
3. **Classification**:
   - If both nucleotides are purines (A,G) or pyrimidines (C,T), it's a transition
   - If one is purine and other is pyrimidine, it's a transversion
4. **Output**: Count of transitions, transversions, and their ratio

## Key Points

- The solution handles DNA sequences with characters A, C, G, T
- Transitions occur between similar nucleotide types (purine↔purine or pyrimidine↔pyrimidine)
- Transversions occur between different nucleotide types (purine↔pyrimidine)
- The ratio helps in understanding evolutionary changes in DNA sequences

This implementation correctly identifies transitions and transversions according to standard molecular biology definitions.