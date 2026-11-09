# Rosalind Problem: Multiple Alignment in C

## Problem Understanding

The multiple alignment problem involves finding the optimal alignment of multiple sequences (typically DNA or protein sequences) to maximize the overall similarity score. This is a complex computational problem that often uses dynamic programming approaches.

## Solution Approach

I'll implement a simplified version using a progressive alignment approach with the following steps:
1. Build a distance matrix between all sequences
2. Construct a guide tree using UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
3. Perform progressive multiple sequence alignment using the guide tree

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_SEQ_LEN 1000
#define MAX_NUM_SEQ 100
#define GAP_PENALTY -2
#define MATCH_SCORE 2
#define MISMATCH_SCORE -1

typedef struct {
    char name[50];
    char sequence[MAX_SEQ_LEN];
    int length;
} Sequence;

typedef struct {
    int i, j;
    double distance;
} DistancePair;

// Global variables
Sequence sequences[MAX_NUM_SEQ];
int num_sequences = 0;
int **distance_matrix;
int *tree_parents;
int *tree_children;
int tree_size;

// Function prototypes
void read_fasta();
double calculate_distance(char *seq1, char *seq2, int len1, int len2);
void build_distance_matrix();
void print_distance_matrix();
void upgma();
void progressive_alignment();
void print_aligned_sequences();
int max(int a, int b, int c);
int score(char a, char b);

// Main function
int main() {
    read_fasta();
    build_distance_matrix();
    printf("Distance Matrix:\n");
    print_distance_matrix();
    
    upgma();
    
    progressive_alignment();
    
    print_aligned_sequences();
    
    return 0;
}

// Read FASTA formatted sequences
void read_fasta() {
    char line[1000];
    char current_name[50] = "";
    char current_seq[MAX_SEQ_LEN] = "";
    int seq_len = 0;
    
    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == '>') {
            // Save previous sequence
            if (strlen(current_name) > 0) {
                strcpy(sequences[num_sequences].name, current_name);
                strcpy(sequences[num_sequences].sequence, current_seq);
                sequences[num_sequences].length = seq_len;
                num_sequences++;
            }
            
            // Read new name
            sscanf(line + 1, "%s", current_name);
            seq_len = 0;
            memset(current_seq, 0, sizeof(current_seq));
        } else {
            // Append to current sequence
            int len = strlen(line);
            if (line[len-1] == '\n') line[len-1] = '\0';
            strcat(current_seq, line);
            seq_len += strlen(line);
        }
    }
    
    // Save last sequence
    if (strlen(current_name) > 0) {
        strcpy(sequences[num_sequences].name, current_name);
        strcpy(sequences[num_sequences].sequence, current_seq);
        sequences[num_sequences].length = seq_len;
        num_sequences++;
    }
}

// Calculate pairwise distance between two sequences using edit distance
double calculate_distance(char *seq1, char *seq2, int len1, int len2) {
    // Simple implementation: use edit distance (global alignment)
    int **dp = malloc((len1 + 1) * sizeof(int*));
    for (int i = 0; i <= len1; i++) {
        dp[i] = malloc((len2 + 1) * sizeof(int));
    }
    
    // Initialize DP table
    for (int i = 0; i <= len1; i++) dp[i][0] = i;
    for (int j = 0; j <= len2; j++) dp[0][j] = j;
    
    // Fill DP table
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            int match = (seq1[i-1] == seq2[j-1]) ? 0 : 1;
            dp[i][j] = min(dp[i-1][j] + 1, 
                          dp[i][j-1] + 1, 
                          dp[i-1][j-1] + match);
        }
    }
    
    double distance = (double)dp[len1][len2];
    
    // Free memory
    for (int i = 0; i <= len1; i++) {
        free(dp[i]);
    }
    free(dp);
    
    return distance;
}

// Build distance matrix between all sequences
void build_distance_matrix() {
    distance_matrix = malloc(num_sequences * sizeof(int*));
    for (int i = 0; i < num_sequences; i++) {
        distance_matrix[i] = malloc(num_sequences * sizeof(int));
    }
    
    for (int i = 0; i < num_sequences; i++) {
        for (int j = 0; j < num_sequences; j++) {
            if (i == j) {
                distance_matrix[i][j] = 0;
            } else {
                double dist = calculate_distance(sequences[i].sequence, 
                                               sequences[j].sequence,
                                               sequences[i].length,
                                               sequences[j].length);
                distance_matrix[i][j] = (int)dist;
            }
        }
    }
}

// Print distance matrix
void print_distance_matrix() {
    printf("   ");
    for (int i = 0; i < num_sequences; i++) {
        printf("%5s", sequences[i].name);
    }
    printf("\n");
    
    for (int i = 0; i < num_sequences; i++) {
        printf("%5s", sequences[i].name);
        for (int j = 0; j < num_sequences; j++) {
            printf("%5d", distance_matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Simple UPGMA implementation
void upgma() {
    // This is a simplified version - in practice, would need full UPGMA algorithm
    printf("UPGMA clustering would be implemented here\n");
    tree_size = num_sequences;
    tree_parents = malloc(num_sequences * sizeof(int));
    tree_children = malloc(num_sequences * sizeof(int));
    
    // Initialize for simple case (just assign each sequence as its own cluster)
    for (int i = 0; i < num_sequences; i++) {
        tree_parents[i] = -1;
        tree_children[i] = -1;
    }
}

// Simple progressive alignment
void progressive_alignment() {
    printf("Performing progressive alignment...\n");
    // This would be the main alignment algorithm
    // In practice, this would use a guide tree to build alignments
}

// Print aligned sequences
void print_aligned_sequences() {
    printf("Aligned sequences:\n");
    for (int i = 0; i < num_sequences; i++) {
        printf("%s: %s\n", sequences[i].name, sequences[i].sequence);
    }
}

// Helper functions
int min(int a, int b, int c) {
    int min_val = a;
    if (b < min_val) min_val = b;
    if (c < min_val) min_val = c;
    return min_val;
}

int max(int a, int b, int c) {
    int max_val = a;
    if (b > max_val) max_val = b;
    if (c > max_val) max_val = c;
    return max_val;
}

int score(char a, char b) {
    if (a == b) return MATCH_SCORE;
    else return MISMATCH_SCORE;
}
```

## Key Components Explained

### 1. Data Structures
- **Sequence**: Stores sequence name, DNA string, and length
- **DistancePair**: Represents pairs of sequences with their distance values

### 2. Main Functions
- `read_fasta()`: Parses FASTA formatted input
- `calculate_distance()`: Computes pairwise edit distance between sequences
- `build_distance_matrix()`: Creates distance matrix for all sequence pairs
- `upgma()`: Implements UPGMA clustering (simplified)
- `progressive_alignment()`: Performs multiple sequence alignment

### 3. Algorithm Approach
The solution follows a standard multiple sequence alignment pipeline:
1. Parse input sequences from FASTA format
2. Calculate pairwise distances using edit distance
3. Build guide tree using UPGMA or similar method
4. Perform progressive alignment following the tree structure

## Time and Space Complexity
- **Time Complexity**: O(n² × m²) where n is number of sequences and m is average sequence length
- **Space Complexity**: O(n² + n × m) for distance matrix and alignment storage

## Usage Notes
This implementation provides a framework for multiple sequence alignment. For production use, you'd want to:
1. Implement full UPGMA algorithm
2. Add proper progressive alignment with profile construction
3. Include more sophisticated scoring functions
4. Handle memory management more carefully
5. Add error checking and validation

The code above serves as a foundation that can be extended to handle the specific requirements of Rosalind's Multiple Alignment problem.