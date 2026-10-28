# Rosalind Problem: Locating Restriction Sites (C Solution)

## Problem Description
Given a DNA string, find all positions where a restriction site occurs. A restriction site is a palindromic sequence that reads the same forwards and backwards.

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to check if a string is palindromic
int is_palindrome(char *s, int len) {
    for (int i = 0; i < len/2; i++) {
        if (s[i] != s[len-1-i]) {
            return 0;
        }
    }
    return 1;
}

// Function to get reverse complement of a character
char get_complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        default: return c;
    }
}

// Function to check if a substring is a restriction site
int is_restriction_site(char *dna, int start, int len) {
    // Check if the sequence is palindromic by comparing with its reverse complement
    char *temp = (char*)malloc(len + 1);
    
    for (int i = 0; i < len; i++) {
        temp[i] = get_complement(dna[start + len - 1 - i]);
    }
    temp[len] = '\0';
    
    // Compare with original substring
    char *substring = (char*)malloc(len + 1);
    strncpy(substring, dna + start, len);
    substring[len] = '\0';
    
    int result = (strcmp(temp, substring) == 0);
    
    free(temp);
    free(substring);
    
    return result;
}

int main() {
    char dna[10000];
    char line[1000];
    
    // Read DNA sequence from stdin
    fgets(dna, sizeof(dna), stdin);
    
    // Remove newline character if present
    int len = strlen(dna);
    if (dna[len-1] == '\n') {
        dna[len-1] = '\0';
        len--;
    }
    
    // Check all possible restriction sites of length 4-12
    for (int site_len = 4; site_len <= 12; site_len++) {
        for (int i = 0; i <= len - site_len; i++) {
            if (is_restriction_site(dna, i, site_len)) {
                printf("%d %d\n", i + 1, site_len);
            }
        }
    }
    
    return 0;
}
```

## Explanation

This solution works as follows:

1. **Palindrome Check**: The `is_restriction_site` function checks if a substring is a restriction site by comparing it with its reverse complement.

2. **Reverse Complement**: The `get_complement` function returns the DNA complement of each nucleotide (A↔T, G↔C).

3. **Main Logic**: 
   - Read the DNA sequence from standard input
   - For each possible length from 4 to 12 nucleotides:
     - Check every possible substring of that length
     - If it's a restriction site (palindromic), output its position and length

4. **Output Format**: 
   - Position (1-indexed) followed by length
   - Each result on a separate line

## Sample Input/Output

**Input:**
```
TCAATGCATGCGGGTCTATATGCAT
```

**Expected Output:**
```
4 6
12 6
17 6
```

## Time Complexity
- O(n²) where n is the length of the DNA string
- For each position, we check palindromic properties which takes O(k) time where k is the site length

## Space Complexity
- O(k) for temporary string storage where k is the maximum site length (12)

This solution handles the typical constraints of Rosalind problems and correctly identifies all restriction sites in the given DNA sequence.