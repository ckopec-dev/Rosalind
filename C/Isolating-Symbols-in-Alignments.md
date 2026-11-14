# Rosalind Problem: Isolating Symbols in Alignments (C Solution)

## Problem Understanding

The problem asks us to find all unique symbols that appear in a multiple sequence alignment and count how many times each symbol appears across all sequences.

## Approach

1. Read the input FASTA format sequences
2. Parse the sequences to extract symbols
3. Count occurrences of each symbol
4. Output the results

## C Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_SEQ_LEN 10000
#define MAX_NUM_SEQS 1000

// Structure to store symbol counts
typedef struct {
    char symbol;
    int count;
} SymbolCount;

// Function to compare symbols for sorting
int compare_symbols(const void *a, const void *b) {
    return ((SymbolCount*)a)->symbol - ((SymbolCount*)b)->symbol;
}

int main() {
    char sequences[MAX_NUM_SEQS][MAX_SEQ_LEN];
    int num_sequences = 0;
    int max_length = 0;
    
    // Read all sequences
    char line[MAX_SEQ_LEN];
    char current_seq[MAX_SEQ_LEN] = "";
    int current_len = 0;
    
    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == '>') {
            // New sequence header
            if (current_len > 0) {
                sequences[num_sequences][current_len] = '\0';
                num_sequences++;
                current_len = 0;
            }
            continue;
        } else {
            // Sequence data
            int len = strlen(line);
            if (len > 0 && line[len-1] == '\n') {
                line[len-1] = '\0';
                len--;
            }
            
            // Append to current sequence
            strcat(current_seq, line);
            current_len += len;
        }
    }
    
    // Don't forget the last sequence
    if (current_len > 0) {
        sequences[num_sequences][current_len] = '\0';
        num_sequences++;
    }
    
    // Find maximum sequence length
    for (int i = 0; i < num_sequences; i++) {
        int len = strlen(sequences[i]);
        if (len > max_length) {
            max_length = len;
        }
    }
    
    // Count symbols at each position
    SymbolCount *symbol_counts = malloc(max_length * sizeof(SymbolCount));
    int symbol_count_size = 0;
    
    // Initialize arrays to track which symbols are seen at each position
    char **seen_symbols = malloc(max_length * sizeof(char*));
    for (int i = 0; i < max_length; i++) {
        seen_symbols[i] = malloc(256 * sizeof(char)); // ASCII characters
        memset(seen_symbols[i], 0, 256);
    }
    
    // Count symbols at each position
    for (int pos = 0; pos < max_length; pos++) {
        int counts[256] = {0};
        
        for (int seq = 0; seq < num_sequences; seq++) {
            if (pos < strlen(sequences[seq])) {
                char symbol = sequences[seq][pos];
                counts[(unsigned char)symbol]++;
                
                // Mark that we've seen this symbol at this position
                if (!seen_symbols[pos][(unsigned char)symbol]) {
                    seen_symbols[pos][(unsigned char)symbol] = 1;
                    symbol_counts[symbol_count_size].symbol = symbol;
                    symbol_counts[symbol_count_size].count = 0;
                    symbol_count_size++;
                }
            }
        }
        
        // Update counts for each symbol at this position
        for (int i = 0; i < symbol_count_size; i++) {
            if (counts[(unsigned char)symbol_counts[i].symbol] > 0) {
                symbol_counts[i].count += counts[(unsigned char)symbol_counts[i].symbol];
            }
        }
    }
    
    // Sort symbols
    qsort(symbol_counts, symbol_count_size, sizeof(SymbolCount), compare_symbols);
    
    // Print results
    for (int i = 0; i < symbol_count_size; i++) {
        printf("%c %d\n", symbol_counts[i].symbol, symbol_counts[i].count);
    }
    
    // Free memory
    for (int i = 0; i < max_length; i++) {
        free(seen_symbols[i]);
    }
    free(seen_symbols);
    free(symbol_counts);
    
    return 0;
}
```

## Alternative Simpler Solution

If we interpret the problem as simply counting all symbols that appear anywhere in the alignment:

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ_LEN 10000
#define MAX_NUM_SEQS 1000

int main() {
    char sequences[MAX_NUM_SEQS][MAX_SEQ_LEN];
    int num_sequences = 0;
    
    // Read sequences
    char line[MAX_SEQ_LEN];
    char current_seq[MAX_SEQ_LEN] = "";
    int current_len = 0;
    
    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == '>') {
            if (current_len > 0) {
                sequences[num_sequences][current_len] = '\0';
                num_sequences++;
                current_len = 0;
            }
        } else {
            int len = strlen(line);
            if (len > 0 && line[len-1] == '\n') {
                line[len-1] = '\0';
                len--;
            }
            strcat(current_seq, line);
            current_len += len;
        }
    }
    
    // Don't forget the last sequence
    if (current_len > 0) {
        sequences[num_sequences][current_len] = '\0';
        num_sequences++;
    }
    
    // Count all symbols
    int symbol_count[256] = {0};
    
    for (int i = 0; i < num_sequences; i++) {
        int len = strlen(sequences[i]);
        for (int j = 0; j < len; j++) {
            symbol_count[(unsigned char)sequences[i][j]]++;
        }
    }
    
    // Print symbols that appear at least once
    for (int i = 0; i < 256; i++) {
        if (symbol_count[i] > 0) {
            printf("%c %d\n", (char)i, symbol_count[i]);
        }
    }
    
    return 0;
}
```

## How to Compile and Run

```bash
gcc -o isolating_symbols isolating_symbols.c
./isolating_symbols < input.txt
```

## Input Format

The input is in FASTA format:
```
>seq1
ACGTACGT
>seq2
ACGTACGT
>seq3
ACGTACGT
```

## Output Format

Symbols sorted alphabetically with their counts:
```
A 6
C 6
G 6
T 6
```

The second solution is simpler and more direct for the typical interpretation of this Rosalind problem, where we want to count all occurrences of each symbol across all sequences in the alignment.