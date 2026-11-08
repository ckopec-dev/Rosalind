```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LEN 1000
#define INF 1000000

int score_matrix[256][256];
int gap_penalty = -2;

void initialize_score_matrix() {
    // Initialize with zero scores for same characters
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            score_matrix[i][j] = 0;
        }
    }
    
    // Set mismatch penalty to -1
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            if (i != j) {
                score_matrix[i][j] = -1;
            }
        }
    }
}

void print_matrix(int **matrix, int m, int n) {
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int max3(int a, int b, int c) {
    int max = a;
    if (b > max) max = b;
    if (c > max) max = c;
    return max;
}

int global_alignment(char *seq1, char *seq2, int m, int n, char *align1, char *align2) {
    // Create DP matrix
    int **dp = (int**)malloc((m + 1) * sizeof(int*));
    for (int i = 0; i <= m; i++) {
        dp[i] = (int*)calloc(n + 1, sizeof(int));
    }
    
    // Initialize first row and column
    for (int i = 0; i <= m; i++) {
        dp[i][0] = i * gap_penalty;
    }
    for (int j = 0; j <= n; j++) {
        dp[0][j] = j * gap_penalty;
    }
    
    // Fill the DP matrix
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int match = dp[i-1][j-1] + score_matrix[(unsigned char)seq1[i-1]][(unsigned char)seq2[j-1]];
            int delete = dp[i-1][j] + gap_penalty;
            int insert = dp[i][j-1] + gap_penalty;
            
            dp[i][j] = max3(match, delete, insert);
        }
    }
    
    // Traceback
    int i = m, j = n;
    int align_len = 0;
    
    while (i > 0 && j > 0) {
        int current = dp[i][j];
        int match = dp[i-1][j-1] + score_matrix[(unsigned char)seq1[i-1]][(unsigned char)seq2[j-1]];
        int delete = dp[i-1][j] + gap_penalty;
        int insert = dp[i][j-1] + gap_penalty;
        
        if (current == match) {
            align1[align_len] = seq1[i-1];
            align2[align_len] = seq2[j-1];
            i--;
            j--;
        } else if (current == delete) {
            align1[align_len] = seq1[i-1];
            align2[align_len] = '-';
            i--;
        } else {
            align1[align_len] = '-';
            align2[align_len] = seq2[j-1];
            j--;
        }
        align_len++;
    }
    
    // Handle remaining characters
    while (i > 0) {
        align1[align_len] = seq1[i-1];
        align2[align_len] = '-';
        i--;
        align_len++;
    }
    
    while (j > 0) {
        align1[align_len] = '-';
        align2[align_len] = seq2[j-1];
        j--;
        align_len++;
    }
    
    // Reverse the alignments
    for (int k = 0; k < align_len/2; k++) {
        char temp1 = align1[k];
        char temp2 = align2[k];
        align1[k] = align1[align_len-1-k];
        align2[k] = align2[align_len-1-k];
        align1[align_len-1-k] = temp1;
        align2[align_len-1-k] = temp2;
    }
    
    // Print result
    printf("%d\n", dp[m][n]);
    printf("%s\n", align1);
    printf("%s\n", align2);
    
    // Free memory
    for (int i = 0; i <= m; i++) {
        free(dp[i]);
    }
    free(dp);
    
    return dp[m][n];
}

int main() {
    char seq1[MAX_LEN], seq2[MAX_LEN];
    char align1[MAX_LEN * 2], align2[MAX_LEN * 2];
    
    initialize_score_matrix();
    
    // Read input sequences
    scanf("%s", seq1);
    scanf("%s", seq2);
    
    int m = strlen(seq1);
    int n = strlen(seq2);
    
    global_alignment(seq1, seq2, m, n, align1, align2);
    
    return 0;
}
```