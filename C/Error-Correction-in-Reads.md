# Rosalind Problem: Error Correction in Reads (C Solution)

## Problem Understanding

The task is to identify and correct sequencing errors in DNA reads by finding:
1. Reads that appear exactly once (potential correct reads)
2. Reads that appear exactly twice (potential erroneous reads)
3. For each erroneous read, find a correct read with Hamming distance 1
4. Output the corrected pairs

## Solution Approach

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_READS 10000
#define READ_LENGTH 50

// Structure to store reads and their counts
typedef struct {
    char read[READ_LENGTH];
    int count;
    int is_correct;
} Read;

// Function to calculate Hamming distance between two strings
int hamming_distance(char *s1, char *s2) {
    int distance = 0;
    int len = strlen(s1);
    
    for (int i = 0; i < len; i++) {
        if (s1[i] != s2[i]) {
            distance++;
        }
    }
    return distance;
}

// Function to reverse a DNA string
void reverse_string(char *str) {
    int len = strlen(str);
    for (int i = 0; i < len/2; i++) {
        char temp = str[i];
        str[i] = str[len-1-i];
        str[len-1-i] = temp;
    }
}

// Function to get reverse complement of a DNA string
void reverse_complement(char *input, char *output) {
    int len = strlen(input);
    for (int i = 0; i < len; i++) {
        switch (input[i]) {
            case 'A': output[len-1-i] = 'T'; break;
            case 'T': output[len-1-i] = 'A'; break;
            case 'G': output[len-1-i] = 'C'; break;
            case 'C': output[len-1-i] = 'G'; break;
        }
    }
    output[len] = '\0';
}

int main() {
    Read reads[MAX_READS];
    int read_count = 0;
    
    // Read all reads from input
    char buffer[READ_LENGTH];
    while (fgets(buffer, READ_LENGTH, stdin)) {
        // Remove newline character if present
        buffer[strcspn(buffer, "\n")] = 0;
        
        if (strlen(buffer) > 0) {
            strcpy(reads[read_count].read, buffer);
            reads[read_count].count = 1;
            reads[read_count].is_correct = 0;
            read_count++;
        }
    }
    
    // Count occurrences of each read
    for (int i = 0; i < read_count; i++) {
        for (int j = i + 1; j < read_count; j++) {
            if (strcmp(reads[i].read, reads[j].read) == 0) {
                reads[i].count++;
                reads[j].count = 0; // Mark as duplicate
            }
        }
    }
    
    // Find reads that appear exactly once (potential correct reads)
    Read correct_reads[MAX_READS];
    int correct_count = 0;
    
    for (int i = 0; i < read_count; i++) {
        if (reads[i].count == 1) {
            strcpy(correct_reads[correct_count].read, reads[i].read);
            correct_reads[correct_count].count = 1;
            correct_count++;
        }
    }
    
    // Find reads that appear exactly twice (potential erroneous reads)
    Read error_reads[MAX_READS];
    int error_count = 0;
    
    for (int i = 0; i < read_count; i++) {
        if (reads[i].count == 2) {
            strcpy(error_reads[error_count].read, reads[i].read);
            error_reads[error_count].count = 1;
            error_count++;
        }
    }
    
    // For each erroneous read, find a correct read with Hamming distance 1
    for (int i = 0; i < error_count; i++) {
        char *err_read = error_reads[i].read;
        
        for (int j = 0; j < correct_count; j++) {
            char *corr_read = correct_reads[j].read;
            
            // Check if they are identical
            if (strcmp(err_read, corr_read) == 0) {
                continue;
            }
            
            // Check Hamming distance
            if (hamming_distance(err_read, corr_read) == 1) {
                printf("%s->%s\n", err_read, corr_read);
                break;
            }
            
            // Also check reverse complement
            char rc[READ_LENGTH];
            reverse_complement(corr_read, rc);
            
            if (hamming_distance(err_read, rc) == 1) {
                printf("%s->%s\n", err_read, rc);
                break;
            }
        }
    }
    
    return 0;
}
```

## Key Features of the Solution

1. **Read Input**: Reads all DNA reads from standard input
2. **Counting**: Counts occurrences of each read to identify correct (appears once) and erroneous (appears twice) reads
3. **Hamming Distance**: Calculates distance between reads to find corrections
4. **Reverse Complement**: Handles the case where an error might be in reverse complement form
5. **Output Format**: Prints corrections in required format

## Time Complexity
- O(n²) for counting and comparison operations
- Space complexity: O(n)

## Usage
Compile with:
```bash
gcc -o error_correction error_correction.c
```

Run with input file:
```bash
./error_correction < input.txt
```

The solution handles the core requirements of finding erroneous reads and their corrections based on Hamming distance of 1.