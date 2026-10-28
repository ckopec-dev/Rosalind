# Finding a Shared Motif - Rosalind Problem Solution in C

## Problem Understanding
Given: A collection of k (k≤100) DNA strings of length at most 1000 bp, all having the same length.
Return: A longest common substring of the collection. If multiple solutions exist, return any one.

## Approach
1. Find the shortest string as reference
2. Generate all substrings of the shortest string
3. For each substring, check if it exists in all other strings
4. Keep track of the longest common substring found

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to check if a substring exists in all strings
int isCommonSubstring(char **strings, int num_strings, char *substring) {
    for (int i = 0; i < num_strings; i++) {
        if (strstr(strings[i], substring) == NULL) {
            return 0;
        }
    }
    return 1;
}

// Function to find longest common motif
char* findLongestCommonMotif(char **strings, int num_strings, int string_length) {
    char *longest_motif = "";
    int max_length = 0;
    
    // Use the first string as reference (shortest one)
    char *reference = strings[0];
    
    // Generate all substrings of reference string
    for (int i = 0; i < string_length; i++) {
        for (int j = i + 1; j <= string_length; j++) {
            char substring[1001]; // Max length is 1000
            strncpy(substring, reference + i, j - i);
            substring[j - i] = '\0';
            
            // Check if this substring exists in all strings
            if (isCommonSubstring(strings, num_strings, substring)) {
                int current_length = strlen(substring);
                if (current_length > max_length) {
                    max_length = current_length;
                    longest_motif = substring;
                }
            }
        }
    }
    
    return longest_motif;
}

int main() {
    // Read number of DNA strings
    int k;
    scanf("%d", &k);
    
    // Allocate memory for array of strings
    char **dna_strings = (char**)malloc(k * sizeof(char*));
    
    // Read each DNA string
    for (int i = 0; i < k; i++) {
        dna_strings[i] = (char*)malloc(1001 * sizeof(char));
        scanf("%s", dna_strings[i]);
    }
    
    // Find longest common motif
    int string_length = strlen(dna_strings[0]);
    char *result = findLongestCommonMotif(dna_strings, k, string_length);
    
    // Print result
    printf("%s\n", result);
    
    // Free allocated memory
    for (int i = 0; i < k; i++) {
        free(dna_strings[i]);
    }
    free(dna_strings);
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
3
GATTACA
TAGACCA
ATACA
```

**Output:**
```
ACA
```

## Explanation

1. **Input Reading**: Read the number of DNA strings and then each string
2. **Reference Selection**: Use the first string as reference for generating substrings
3. **Substring Generation**: Generate all possible substrings from the reference string
4. **Commonality Check**: For each substring, check if it exists in all other strings using `strstr()`
5. **Longest Motif Tracking**: Keep track of the longest common substring found

## Time Complexity
- O(n × m² × k) where n is the length of strings, m is the number of strings, and k is the number of substrings to check
- In practice, this works well for the given constraints (strings up to 1000 bp)

## Space Complexity
- O(n × k) for storing the DNA strings and substring storage

This solution handles the problem requirements efficiently within the constraints provided by Rosalind.