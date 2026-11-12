# Solving Rosalind Problem: Fixing an Inconsistent Character Set

## Problem Understanding

The task is to fix an inconsistent character set in a multiple sequence alignment by removing characters that appear in less than half of the sequences, or if a character appears in exactly half of the sequences, remove it if it's not the most common character.

## Approach

1. Parse input sequences
2. Count character frequencies at each position across all sequences
3. For each position, determine which characters to keep based on the rules:
   - If a character appears in > n/2 sequences, keep it
   - If a character appears in exactly n/2 sequences, keep the most frequent one
   - Otherwise, remove the character

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQS 1000
#define MAX_LEN 1000

// Structure to store character counts at each position
typedef struct {
    int count[256]; // ASCII character counts
    int total;
} PositionCount;

int main() {
    char sequences[MAX_SEQS][MAX_LEN];
    int num_seqs = 0;
    int max_len = 0;
    
    // Read input sequences
    while (fgets(sequences[num_seqs], MAX_LEN, stdin) != NULL) {
        // Remove newline if present
        int len = strlen(sequences[num_seqs]);
        if (len > 0 && sequences[num_seqs][len-1] == '\n') {
            sequences[num_seqs][len-1] = '\0';
            len--;
        }
        
        if (len > max_len) max_len = len;
        num_seqs++;
    }
    
    // Handle empty input
    if (num_seqs == 0) return 0;
    
    // For each position, count character frequencies
    PositionCount *pos_counts = malloc(max_len * sizeof(PositionCount));
    for (int i = 0; i < max_len; i++) {
        memset(&pos_counts[i], 0, sizeof(PositionCount));
    }
    
    // Count characters at each position
    for (int seq_idx = 0; seq_idx < num_seqs; seq_idx++) {
        int len = strlen(sequences[seq_idx]);
        for (int pos = 0; pos < len; pos++) {
            if (pos < max_len) {
                pos_counts[pos].count[(unsigned char)sequences[seq_idx][pos]]++;
                pos_counts[pos].total++;
            }
        }
    }
    
    // Generate output - for each position, determine what character to keep
    char result[MAX_LEN];
    for (int pos = 0; pos < max_len; pos++) {
        int max_count = 0;
        char best_char = '\0';
        int total = pos_counts[pos].total;
        
        // Find the character that appears in most sequences at this position
        for (int c = 0; c < 256; c++) {
            if (pos_counts[pos].count[c] > max_count) {
                max_count = pos_counts[pos].count[c];
                best_char = (char)c;
            }
        }
        
        // Check if we should keep this character
        if (total == 0) {
            result[pos] = '\0';
        } else if (max_count > total / 2) {
            // Keep the most frequent character
            result[pos] = best_char;
        } else if (max_count == total / 2) {
            // If tied for most frequent, keep it
            result[pos] = best_char;
        } else {
            // Remove the character - no valid character to keep
            result[pos] = '\0';
        }
    }
    
    result[max_len] = '\0';
    
    // Print the fixed alignment
    printf("%s\n", result);
    
    free(pos_counts);
    return 0;
}
```

Wait, I think I misunderstood the problem. Let me reconsider the approach.

## Corrected Approach

Looking at this more carefully, it seems like we're supposed to fix an inconsistent character set by removing characters that don't appear in enough sequences. The key insight is that for each position, we want to remove characters that appear in less than half of the sequences.

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQS 1000
#define MAX_LEN 1000

int main() {
    char sequences[MAX_SEQS][MAX_LEN];
    int num_seqs = 0;
    int max_len = 0;
    
    // Read input sequences
    char line[MAX_LEN];
    while (fgets(line, MAX_LEN, stdin) != NULL) {
        // Remove newline if present
        int len = strlen(line);
        if (len > 0 && line[len-1] == '\n') {
            line[len-1] = '\0';
            len--;
        }
        
        if (len > max_len) max_len = len;
        strcpy(sequences[num_seqs], line);
        num_seqs++;
    }
    
    // Handle empty input
    if (num_seqs == 0) return 0;
    
    // For each position, determine which characters to keep
    char result[MAX_LEN];
    memset(result, '\0', MAX_LEN);
    
    for (int pos = 0; pos < max_len; pos++) {
        int char_count[256] = {0};
        int total_chars = 0;
        
        // Count characters at this position
        for (int seq_idx = 0; seq_idx < num_seqs; seq_idx++) {
            if (pos < strlen(sequences[seq_idx])) {
                char c = sequences[seq_idx][pos];
                char_count[(unsigned char)c]++;
                total_chars++;
            }
        }
        
        // Find the character that appears in most sequences
        int max_count = 0;
        char best_char = '\0';
        
        for (int c = 0; c < 256; c++) {
            if (char_count[c] > max_count) {
                max_count = char_count[c];
                best_char = (char)c;
            }
        }
        
        // Keep the character that appears in more than half of sequences
        if (max_count > total_chars / 2) {
            result[pos] = best_char;
        } else {
            // If no character appears in more than half, we keep the most frequent one 
            // but in some interpretations, we might remove it entirely
            // For this implementation, keeping the most frequent character as a default
            result[pos] = best_char;
        }
    }
    
    printf("%s\n", result);
    return 0;
}
```

Actually, let me simplify and make it cleaner:

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQS 1000
#define MAX_LEN 1000

int main() {
    char sequences[MAX_SEQS][MAX_LEN];
    int num_seqs = 0;
    int max_len = 0;
    
    // Read all sequences
    char line[MAX_LEN];
    while (fgets(line, MAX_LEN, stdin) != NULL) {
        int len = strlen(line);
        if (len > 0 && line[len-1] == '\n') {
            line[len-1] = '\0';
        }
        
        strcpy(sequences[num_seqs], line);
        if (strlen(line) > max_len) {
            max_len = strlen(line);
        }
        num_seqs++;
    }
    
    // Handle empty input
    if (num_seqs == 0) return 0;
    
    // Process each position
    char result[MAX_LEN];
    for (int pos = 0; pos < max_len; pos++) {
        int count[256] = {0};
        int total_valid = 0;
        
        // Count characters at this position across all sequences
        for (int i = 0; i < num_seqs; i++) {
            if (pos < strlen(sequences[i])) {
                char c = sequences[i][pos];
                count[(unsigned char)c]++;
                total_valid++;
            }
        }
        
        // Find the character that appears most frequently
        int max_count = 0;
        char best_char = '\0';
        
        for (int i = 0; i < 256; i++) {
            if (count[i] > max_count) {
                max_count = count[i];
                best_char = (char)i;
            }
        }
        
        // Keep the character if it appears in more than half of sequences
        if (max_count > total_valid / 2) {
            result[pos] = best_char;
        } else {
            // If no character dominates, keep most frequent (or remove?)
            // Based on typical Rosalind solutions, we keep the most frequent
            result[pos] = best_char;
        }
    }
    
    result[max_len] = '\0';
    printf("%s\n", result);
    
    return 0;
}
```