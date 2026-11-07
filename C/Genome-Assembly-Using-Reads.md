# Rosalind Problem: Genome Assembly Using Reads (C Solution)

## Problem Understanding

Given a collection of DNA reads, we need to find the shortest possible DNA string that contains all the reads as substrings. This is essentially finding the shortest superstring that contains all given strings.

## Approach

I'll use a greedy approach:
1. Find the maximum overlap between any two reads
2. Merge the pair with maximum overlap
3. Repeat until only one string remains

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_READS 1000
#define MAX_LENGTH 1000

// Function to calculate overlap between two strings
int get_overlap(char* a, char* b) {
    int len_a = strlen(a);
    int len_b = strlen(b);
    int max_overlap = 0;
    
    // Try all possible overlaps
    for (int i = 1; i <= len_a && i <= len_b; i++) {
        if (strncmp(a + len_a - i, b, i) == 0) {
            if (i > max_overlap) {
                max_overlap = i;
            }
        }
    }
    
    return max_overlap;
}

// Function to merge two strings with maximum overlap
char* merge_strings(char* a, char* b, int overlap) {
    int len_a = strlen(a);
    int len_b = strlen(b);
    int total_len = len_a + len_b - overlap;
    
    char* result = (char*)malloc((total_len + 1) * sizeof(char));
    
    // Copy first string
    strcpy(result, a);
    
    // Append remaining part of second string
    strcat(result, b + overlap);
    
    return result;
}

// Function to find the shortest superstring
char* find_shortest_superstring(char** reads, int num_reads) {
    if (num_reads == 1) {
        char* result = (char*)malloc((strlen(reads[0]) + 1) * sizeof(char));
        strcpy(result, reads[0]);
        return result;
    }
    
    // Create a copy of reads to work with
    char** working_reads = (char**)malloc(num_reads * sizeof(char*));
    for (int i = 0; i < num_reads; i++) {
        working_reads[i] = (char*)malloc((strlen(reads[i]) + 1) * sizeof(char));
        strcpy(working_reads[i], reads[i]);
    }
    
    int current_num = num_reads;
    
    while (current_num > 1) {
        int max_overlap = -1;
        int best_i = -1, best_j = -1;
        char* merged_string = NULL;
        
        // Find the pair with maximum overlap
        for (int i = 0; i < current_num; i++) {
            for (int j = 0; j < current_num; j++) {
                if (i != j) {
                    int overlap = get_overlap(working_reads[i], working_reads[j]);
                    if (overlap > max_overlap) {
                        max_overlap = overlap;
                        best_i = i;
                        best_j = j;
                    }
                }
            }
        }
        
        // If no overlap found, just merge first two
        if (max_overlap <= 0) {
            merged_string = merge_strings(working_reads[0], working_reads[1], 0);
            working_reads[0] = merged_string;
            
            // Shift remaining elements
            for (int i = 1; i < current_num - 1; i++) {
                free(working_reads[i]);
                working_reads[i] = working_reads[i + 1];
            }
            current_num--;
        } else {
            // Merge the best pair
            merged_string = merge_strings(working_reads[best_i], working_reads[best_j], max_overlap);
            
            // Place merged string at position best_i
            free(working_reads[best_i]);
            working_reads[best_i] = merged_string;
            
            // Shift remaining elements to remove element at best_j
            for (int i = best_j; i < current_num - 1; i++) {
                free(working_reads[i]);
                working_reads[i] = working_reads[i + 1];
            }
            current_num--;
        }
    }
    
    char* result = (char*)malloc((strlen(working_reads[0]) + 1) * sizeof(char));
    strcpy(result, working_reads[0]);
    
    // Free memory
    for (int i = 0; i < num_reads; i++) {
        free(working_reads[i]);
    }
    free(working_reads);
    
    return result;
}

int main() {
    // Read input
    char** reads = (char**)malloc(MAX_READS * sizeof(char*));
    int num_reads = 0;
    
    char buffer[MAX_LENGTH];
    
    // Read all reads until EOF
    while (fgets(buffer, MAX_LENGTH, stdin)) {
        // Remove newline character if present
        buffer[strcspn(buffer, "\n")] = 0;
        
        if (strlen(buffer) > 0) {
            reads[num_reads] = (char*)malloc((strlen(buffer) + 1) * sizeof(char));
            strcpy(reads[num_reads], buffer);
            num_reads++;
        }
    }
    
    // Find shortest superstring
    char* result = find_shortest_superstring(reads, num_reads);
    
    printf("%s\n", result);
    
    // Free memory
    for (int i = 0; i < num_reads; i++) {
        free(reads[i]);
    }
    free(reads);
    free(result);
    
    return 0;
}
```

## Explanation

1. **get_overlap()**: Calculates the maximum overlap between two strings by checking all possible overlaps from the end of string A to the beginning of string B.

2. **merge_strings()**: Combines two strings by overlapping them maximally and returns the merged result.

3. **find_shortest_superstring()**: Implements a greedy approach:
   - Continuously finds the pair of strings with maximum overlap
   - Merges them into a single string
   - Repeats until only one string remains

4. **Main function**: 
   - Reads all input reads from standard input
   - Calls the assembly function
   - Outputs the final superstring

## Time Complexity
- O(n² × m) where n is the number of reads and m is the average read length
- This is a greedy approximation algorithm, not guaranteed to find the optimal solution but works well in practice for most cases

## Usage
Compile with: `gcc -o genome_assembly genome_assembly.c`
Run with input: `./genome_assembly < input.txt`