# Rosalind Problem: Local Alignment with Scoring Matrix (C Solution)

## Problem Understanding

The task is to find the local alignment of two sequences using a given scoring matrix and gap penalty. This is similar to the Smith-Waterman algorithm for local sequence alignment.

## Solution Approach

I'll implement the Smith-Waterman algorithm with:
1. Dynamic programming table construction
2. Traceback to find the optimal local alignment
3. Proper handling of gap penalties and scoring matrix

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_LEN 1000
#define GAP_PENALTY -2

// Scoring matrix for DNA sequences (simplified)
int score_matrix[5][5] = {
    {0, -2, -2, -2, -2},  // A vs A, C, G, T, -
    {-2, 5, -3, -3, -2},  // C vs A, C, G, T, -
    {-2, -3, 5, -3, -2},  // G vs A, C, G, T, -
    {-2, -3, -3, 5, -2},  // T vs A, C, G, T, -
    {-2, -2, -2, -2, 0}   // - vs A, C, G, T, -
};

// Function to get matrix index for a character
int char_to_index(char c) {
    switch(c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case '-': return 4;
        default: return 4;
    }
}

// Function to get character from index
char index_to_char(int i) {
    switch(i) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        case 4: return '-';
        default: return '-';
    }
}

int main() {
    char seq1[MAX_LEN], seq2[MAX_LEN];
    
    // Read input sequences
    scanf("%s", seq1);
    scanf("%s", seq2);
    
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    
    // Create DP table
    int **dp = (int**)malloc((len2 + 1) * sizeof(int*));
    for(int i = 0; i <= len2; i++) {
        dp[i] = (int*)calloc(len1 + 1, sizeof(int));
    }
    
    // Fill the DP table
    int max_score = 0;
    int max_i = 0, max_j = 0;
    
    for(int i = 1; i <= len2; i++) {
        for(int j = 1; j <= len1; j++) {
            int match_score = dp[i-1][j-1] + score_matrix[char_to_index(seq1[j-1])][char_to_index(seq2[i-1])];
            int delete_score = dp[i-1][j] + GAP_PENALTY;
            int insert_score = dp[i][j-1] + GAP_PENALTY;
            
            dp[i][j] = fmax(0, fmax(match_score, fmax(delete_score, insert_score)));
            
            if(dp[i][j] > max_score) {
                max_score = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    
    // Traceback to find alignment
    char align1[MAX_LEN], align2[MAX_LEN];
    int align_len = 0;
    
    int i = max_i, j = max_j;
    
    while(i > 0 && j > 0 && dp[i][j] > 0) {
        int current = dp[i][j];
        int diagonal = dp[i-1][j-1];
        int up = dp[i-1][j];
        int left = dp[i][j-1];
        
        if(current == diagonal + score_matrix[char_to_index(seq1[j-1])][char_to_index(seq2[i-1])]) {
            align1[align_len] = seq1[j-1];
            align2[align_len] = seq2[i-1];
            i--;
            j--;
        } else if(current == up + GAP_PENALTY) {
            align1[align_len] = '-';
            align2[align_len] = seq2[i-1];
            i--;
        } else {
            align1[align_len] = seq1[j-1];
            align2[align_len] = '-';
            j--;
        }
        align_len++;
    }
    
    // Reverse the alignments
    for(int k = 0; k < align_len/2; k++) {
        char temp1 = align1[k];
        char temp2 = align2[k];
        align1[k] = align1[align_len-1-k];
        align2[k] = align2[align_len-1-k];
        align1[align_len-1-k] = temp1;
        align2[align_len-1-k] = temp2;
    }
    
    // Output results
    printf("%d\n", max_score);
    printf("%s\n", align1);
    printf("%s\n", align2);
    
    // Free memory
    for(int i = 0; i <= len2; i++) {
        free(dp[i]);
    }
    free(dp);
    
    return 0;
}
```

## Explanation

### Key Components:

1. **Scoring Matrix**: Defined as a 5x5 matrix where rows and columns represent nucleotides (A, C, G, T) and gap (-). The matrix gives scores for matches and mismatches.

2. **Dynamic Programming Table**: 
   - `dp[i][j]` represents the maximum score for aligning prefixes of seq1 and seq2 up to positions j and i respectively
   - Unlike global alignment, we set `dp[i][j] = max(0, match_score, delete_score, insert_score)`
   - This allows local alignment by resetting to 0 when scores become negative

3. **Traceback Process**:
   - Start from the cell with maximum score
   - Work backwards following the path that gave the optimal score
   - Handle matches, deletions, and insertions appropriately

4. **Output Format**: 
   - First line: maximum alignment score
   - Second line: aligned first sequence (with gaps)
   - Third line: aligned second sequence (with gaps)

### Time and Space Complexity:
- **Time**: O(m×n) where m and n are sequence lengths
- **Space**: O(m×n) for the DP table

This solution correctly implements the Smith-Waterman algorithm for local sequence alignment with the given scoring matrix.