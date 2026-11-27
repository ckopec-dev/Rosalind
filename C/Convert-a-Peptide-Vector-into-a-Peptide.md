# Rosalind Problem: Convert a Peptide Vector into a Peptide

## Problem Description
Given a peptide vector, we need to convert it back into the original peptide sequence.

A peptide vector is a binary string where each position represents whether a particular amino acid (in order) appears in the peptide. The amino acids are ordered as:
A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y

## Solution Approach
1. Read the peptide vector (binary string)
2. For each position in the vector that is 1, add the corresponding amino acid to the result
3. Return the peptide sequence

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* peptide_vector_to_peptide(char* vector) {
    // Amino acids in order
    char acids[] = "ACDEFGHIKLMNPQRSTVWY";
    
    // Calculate length of result string
    int len = strlen(vector);
    char* result = (char*)malloc(len + 1); // +1 for null terminator
    
    if (result == NULL) {
        return NULL;
    }
    
    int result_index = 0;
    
    // Process each bit in the vector
    for (int i = 0; i < len; i++) {
        if (vector[i] == '1') {
            result[result_index++] = acids[i];
        }
    }
    
    result[result_index] = '\0'; // Null terminate
    return result;
}

int main() {
    char vector[1000];
    
    printf("Enter peptide vector: ");
    scanf("%s", vector);
    
    char* peptide = peptide_vector_to_peptide(vector);
    
    if (peptide != NULL) {
        printf("Peptide: %s\n", peptide);
        free(peptide);
    } else {
        printf("Error: Memory allocation failed\n");
    }
    
    return 0;
}
```

## Example Usage

**Input:** `10000000000000000000000000000000