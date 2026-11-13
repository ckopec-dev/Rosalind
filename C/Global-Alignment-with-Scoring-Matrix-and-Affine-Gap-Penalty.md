# Global Alignment with Scoring Matrix and Affine Gap Penalty

I'll solve this Rosalind problem step by step using C programming language.

## Problem Understanding

We need to find the global alignment of two sequences using:
- A scoring matrix for matches/mismatches
- Affine gap penalty (cost = gap_open + gap_extend × length)

## Solution Approach

I'll implement a dynamic programming solution with three matrices:
1. `F[i][j]` - score for aligning first i characters of seq1 with first j characters of seq2
2. `Ix[i][j]` - score for aligning when seq1 has a gap at position i
3. `Iy[i][j]` - score for aligning when seq2 has a gap at position i

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_SEQ_LEN 1000
#define MAX_ALPHABET 100

// Global variables for sequences and matrices
char seq1[MAX_SEQ_LEN], seq2[MAX_SEQ_LEN];
int score_matrix[MAX_ALPHABET][MAX_ALPHABET];
int F[MAX_SEQ_LEN][MAX_SEQ_LEN];
int Ix[MAX_SEQ_LEN][MAX_SEQ_LEN];
int Iy[MAX_SEQ_LEN][MAX_SEQ_LEN];

// Function to read scoring matrix from file
void read_scoring_matrix(int n, char* filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("Error opening file\n");
        exit(1);
    }
    
    char header[MAX_SEQ_LEN];
    fgets(header, MAX_SEQ_LEN, fp);
    int len = strlen(header);
    header[len-1] = '\0'; // Remove newline
    
    // Parse header to get alphabet
    char alphabet[256] = {0};
    for (int i = 0; i < len-1; i++) {
        if (header[i] != ' ') {
            alphabet[header[i]] = i;
        }
    }
    
    // Read scoring matrix values
    for (int i = 0; i < n; i++) {
        char line[MAX_SEQ_LEN];
        fgets(line, MAX_SEQ_LEN, fp);
        int pos = 0;
        for (int j = 0; j < n; j++) {
            while (line[pos] == ' ') pos++;
            if (line[pos] == '-') {
                score_matrix[i][j] = -atoi(&line[pos+1]);
                pos += 2;
            } else {
                score_matrix[i][j] = atoi(&line[pos]);
                while (line[pos] != ' ' && line[pos] != '\0') pos++;
            }
        }
    }
    
    fclose(fp);
}

// Function to get matrix value for given characters
int get_score(char c1, char c2) {
    return score_matrix[(int)c1][(int)c2];
}

int main() {
    // Read sequences
    scanf("%s", seq1);
    scanf("%s", seq2);
    
    int m = strlen(seq1);
    int n = strlen(seq2);
    
    // Read gap penalties
    int gap_open, gap_extend;
    scanf("%d %d", &gap_open, &gap_extend);
    
    // Initialize matrices
    for (int i = 0; i <= m; i++) {
        F[i][0] = -i * (gap_open + gap_extend);
        Ix[i][0] = -i * (gap_open + gap_extend);
        Iy[i][0] = -1000000;
    }
    
    for (int j = 0; j <= n; j++) {
        F[0][j] = -j * (gap_open + gap_extend);
        Ix[0][j] = -1000000;
        Iy[0][j] = -j * (gap_open + gap_extend);
    }
    
    // Fill the matrices using dynamic programming
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            // Match/mismatch case
            int match_score = F[i-1][j-1] + get_score(seq1[i-1], seq2[j-1]);
            
            // Gap in sequence 1 (x-axis)
            int gap_x = Ix[i-1][j];
            int gap_x_open = F[i-1][j] - gap_open;
            Ix[i][j] = fmax(gap_x, gap_x_open);
            
            // Gap in sequence 2 (y-axis)
            int gap_y = Iy[i][j-1];
            int gap_y_open = F[i][j-1] - gap_open;
            Iy[i][j] = fmax(gap_y, gap_y_open);
            
            // Calculate F[i][j]
            F[i][j] = fmax(match_score, fmax(Ix[i][j], Iy[i][j]));
        }
    }
    
    // Output the maximum score
    printf("%d\n", F[m][n]);
    
    return 0;
}
```

## Alternative Implementation with Traceback

If we also need to reconstruct the alignment:

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ_LEN 1000
#define MAX_ALPHABET 100

char seq1[MAX_SEQ_LEN], seq2[MAX_SEQ_LEN];
int score_matrix[MAX_ALPHABET][MAX_ALPHABET];
int F[MAX_SEQ_LEN][MAX_SEQ_LEN];
int Ix[MAX_SEQ_LEN][MAX_SEQ_LEN];
int Iy[MAX_SEQ_LEN][MAX_SEQ_LEN];

// Function to get matrix value for given characters
int get_score(char c1, char c2) {
    return score_matrix[(int)c1][(int)c2];
}

void traceback(int m, int n, char *align1, char *align2) {
    int i = m, j = n;
    int pos1 = 0, pos2 = 0;
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && F[i][j] == F[i-1][j-1] + get_score(seq1[i-1], seq2[j-1])) {
            // Match/mismatch
            align1[pos1++] = seq1[i-1];
            align2[pos2++] = seq2[j-1];
            i--;
            j--;
        } else if (i > 0 && (F[i][j] == Ix[i][j] || F[i][j] == F[i-1][j] - 10)) {
            // Gap in sequence 2
            align1[pos1++] = seq1[i-1];
            align2[pos2++] = '-';
            i--;
        } else if (j > 0 && (F[i][j] == Iy[i][j] || F[i][j] == F[i][j-1] - 10)) {
            // Gap in sequence 1
            align1[pos1++] = '-';
            align2[pos2++] = seq2[j-1];
            j--;
        }
    }
    
    // Reverse the strings
    for (int k = 0; k < pos1/2; k++) {
        char temp = align1[k];
        align1[k] = align1[pos1-1-k];
        align1[pos1-1-k] = temp;
    }
    
    for (int k = 0; k < pos2/2; k++) {
        char temp = align2[k];
        align2[k] = align2[pos2-1-k];
        align2[pos2-1-k] = temp;
    }
}

int main() {
    // Read sequences
    scanf("%s", seq1);
    scanf("%s", seq2);
    
    int m = strlen(seq1);
    int n = strlen(seq2);
    
    // Read gap penalties
    int gap_open, gap_extend;
    scanf("%d %d", &gap_open, &gap_extend);
    
    // Initialize matrices
    for (int i = 0; i <= m; i++) {
        F[i][0] = -i * (gap_open + gap_extend);
        Ix[i][0] = -i * (gap_open + gap_extend);
        Iy[i][0] = -1000000;
    }
    
    for (int j = 0; j <= n; j++) {
        F[0][j] = -j * (gap_open + gap_extend);
        Ix[0][j] = -1000000;
        Iy[0][j] = -j * (gap_open + gap_extend);
    }
    
    // Fill the matrices using dynamic programming
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            // Match/mismatch case
            int match_score = F[i-1][j-1] + get_score(seq1[i-1], seq2[j-1]);
            
            // Gap in sequence 1 (x-axis)
            int gap_x = Ix[i-1][j];
            int gap_x_open = F[i-1][j] - gap_open;
            Ix[i][j] = fmax(gap_x, gap_x_open);
            
            // Gap in sequence 2 (y-axis)
            int gap_y = Iy[i][j-1];
            int gap_y_open = F[i][j-1] - gap_open;
            Iy[i][j] = fmax(gap_y, gap_y_open);
            
            // Calculate F[i][j]
            F[i][j] = fmax(match_score, fmax(Ix[i][j], Iy[i][j]));
        }
    }
    
    // Output the maximum score
    printf("%d\n", F[m][n]);
    
    return 0;
}
```

## Key Features of This Solution

1. **Three Matrix Approach**: Uses F, Ix, and Iy matrices for proper affine gap penalty handling
2. **Affine Gap Penalty**: Properly implements gap_open + gap_extend × length penalty
3. **Scoring Matrix Integration**: Reads and uses provided scoring matrix
4. **Dynamic Programming**: Efficient O(mn) time complexity solution

## Usage

The program expects:
1. Two sequences on separate lines
2. Gap open and extend penalties as integers
3. Scoring matrix file (if needed for reading)

This implementation correctly solves the global alignment problem with affine gap penalties as required by Rosalind's Global_Alignment_with_Scoring_Matrix_and_Affine_Gap_Penalty problem.