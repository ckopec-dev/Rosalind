# Rosalind Problem: Overlap Graphs (Overlap_Graphs)

## Problem Description
Given a collection of DNA strings, construct a directed graph where each string is a node, and there's a directed edge from string A to string B if the suffix of A (of length k) matches the prefix of B (of length k), where k is typically 3.

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_STRINGS 1000
#define MAX_LENGTH 1000
#define K 3

typedef struct {
    char name[50];
    char sequence[MAX_LENGTH];
} DNA_String;

// Function to check if suffix of string1 matches prefix of string2
int suffix_matches_prefix(char *str1, char *str2, int len) {
    int str1_len = strlen(str1);
    int str2_len = strlen(str2);
    
    if (str1_len < len || str2_len < len) return 0;
    
    // Compare last 'len' characters of str1 with first 'len' characters of str2
    for (int i = 0; i < len; i++) {
        if (str1[str1_len - len + i] != str2[i]) {
            return 0;
        }
    }
    return 1;
}

int main() {
    DNA_String dna_strings[MAX_STRINGS];
    int num_strings = 0;
    
    // Read input
    char line[200];
    char current_name[50];
    char current_seq[MAX_LENGTH] = "";
    
    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == '>') {
            // Save previous sequence if exists
            if (strlen(current_seq) > 0) {
                strcpy(dna_strings[num_strings].name, current_name);
                strcpy(dna_strings[num_strings].sequence, current_seq);
                num_strings++;
            }
            
            // Start new sequence
            line[strlen(line)-1] = '\0'; // Remove newline
            strcpy(current_name, line + 1); // Skip '>'
            strcpy(current_seq, "");
        } else {
            // Append to current sequence
            line[strlen(line)-1] = '\0'; // Remove newline
            strcat(current_seq, line);
        }
    }
    
    // Don't forget the last sequence
    if (strlen(current_seq) > 0) {
        strcpy(dna_strings[num_strings].name, current_name);
        strcpy(dna_strings[num_strings].sequence, current_seq);
        num_strings++;
    }
    
    // Build overlap graph
    for (int i = 0; i < num_strings; i++) {
        for (int j = 0; j < num_strings; j++) {
            if (i != j && suffix_matches_prefix(dna_strings[i].sequence, 
                                               dna_strings[j].sequence, K)) {
                printf("%s %s\n", dna_strings[i].name, dna_strings[j].name);
            }
        }
    }
    
    return 0;
}
```

## Explanation

1. **Data Structure**: We use a `DNA_String` struct to store both the name and sequence of each DNA string.

2. **Input Parsing**: The program reads FASTA formatted input where:
   - Lines starting with '>' indicate new sequences
   - Sequences are built by concatenating all lines between '>' markers

3. **Overlap Detection**: For each pair of DNA strings (i,j), we check if the suffix of string i (last K characters) matches the prefix of string j (first K characters).

4. **Output**: When an overlap is found, we print the names of the two strings in the format "from to".

## Sample Input
```
>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATTTT
>Rosalind_5013
AAATTTT
```

## Sample Output
```
Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_0498 Rosalind_5013
```

## Time Complexity
O(n² × k) where n is the number of DNA strings and k is the overlap length (3).

## Space Complexity
O(n × m) where n is the number of strings and m is the average string length.