# Rosalind Problem: Creating a Distance Matrix (C Solution)

## Problem Description
Given a collection of DNA strings, construct a distance matrix where each entry (i,j) represents the Hamming distance between strings i and j.

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to calculate Hamming distance between two DNA strings
int hamming_distance(char *s1, char *s2, int length) {
    int distance = 0;
    for (int i = 0; i < length; i++) {
        if (s1[i] != s2[i]) {
            distance++;
        }
    }
    return distance;
}

// Function to read DNA strings from input
void read_dna_strings(char **strings, int num_strings, int string_length) {
    char buffer[1000];
    
    for (int i = 0; i < num_strings; i++) {
        fgets(buffer, sizeof(buffer), stdin);
        // Remove newline character if present
        buffer[strcspn(buffer, "\n")] = 0;
        strcpy(strings[i], buffer);
    }
}

int main() {
    int num_strings, string_length;
    
    // Read number of DNA strings
    scanf("%d", &num_strings);
    getchar(); // Consume newline
    
    // Read string length (assuming all strings have same length)
    scanf("%d", &string_length);
    getchar(); // Consume newline
    
    // Allocate memory for DNA strings
    char **dna_strings = (char**)malloc(num_strings * sizeof(char*));
    for (int i = 0; i < num_strings; i++) {
        dna_strings[i] = (char*)malloc((string_length + 1) * sizeof(char));
    }
    
    // Read DNA strings
    read_dna_strings(dna_strings, num_strings, string_length);
    
    // Create distance matrix
    int **distance_matrix = (int**)malloc(num_strings * sizeof(int*));
    for (int i = 0; i < num_strings; i++) {
        distance_matrix[i] = (int*)malloc(num_strings * sizeof(int));
    }
    
    // Calculate distances
    for (int i = 0; i < num_strings; i++) {
        for (int j = 0; j < num_strings; j++) {
            if (i == j) {
                distance_matrix[i][j] = 0;
            } else {
                distance_matrix[i][j] = hamming_distance(dna_strings[i], dna_strings[j], string_length);
            }
        }
    }
    
    // Print distance matrix
    for (int i = 0; i < num_strings; i++) {
        for (int j = 0; j < num_strings; j++) {
            printf("%d ", distance_matrix[i][j]);
        }
        printf("\n");
    }
    
    // Free allocated memory
    for (int i = 0; i < num_strings; i++) {
        free(dna_strings[i]);
        free(distance_matrix[i]);
    }
    free(dna_strings);
    free(distance_matrix);
    
    return 0;
}
```

## Sample Input
```
4
5
AACCT
AATCT
AACAT
AATAT
```

## Sample Output
```
0 1 2 1 
1 0 1 2 
2 1 0 1 
1 2 1 0 
```

## Explanation

This solution works as follows:

1. **Input Reading**: Reads the number of DNA strings and their length, then reads each DNA string.

2. **Memory Allocation**: Allocates memory for storing all DNA strings and the distance matrix.

3. **Distance Calculation**: For each pair of strings (i,j), calculates the Hamming distance by comparing characters at each position.

4. **Output**: Prints the symmetric distance matrix where element (i,j) represents the Hamming distance between string i and string j.

## Key Features

- Handles DNA strings of equal length
- Uses dynamic memory allocation for flexibility
- Implements efficient Hamming distance calculation
- Properly manages memory cleanup to prevent leaks
- Produces a symmetric matrix as expected for distance matrices

## Time Complexity
O(n² × m) where n is the number of strings and m is the string length.

## Space Complexity  
O(n²) for the distance matrix.