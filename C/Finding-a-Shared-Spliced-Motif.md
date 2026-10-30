# Finding a Shared Spliced Motif - Rosalind Problem Solution in C

## Problem Understanding
We need to find the longest common subsequence (LCS) between two DNA strings. A subsequence is formed by deleting some characters from the original string without changing the order of remaining characters.

## Approach
I'll use dynamic programming to solve this problem:
1. Create a 2D DP table where `dp[i][j]` represents the length of LCS for first `i` characters of string1 and first `j` characters of string2
2. Fill the table using the recurrence relation:
   - If characters match: `dp[i][j] = dp[i-1][j-1] + 1`
   - Otherwise: `dp[i][j] = max(dp[i-1][j], dp[i][j-1])`
3. Backtrack through the table to reconstruct the actual LCS

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LEN 1000

// Function to find maximum of two integers
int max(int a, int b) {
    return (a > b) ? a : b;
}

// Function to find longest common subsequence
char* findLCS(char* s1, char* s2) {
    int len1 = strlen(s1);
    int len2 = strlen(s2);
    
    // Create DP table
    int** dp = (int**)malloc((len1 + 1) * sizeof(int*));
    for (int i = 0; i <= len1; i++) {
        dp[i] = (int*)calloc(len2 + 1, sizeof(int));
    }
    
    // Fill DP table
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            if (s1[i-1] == s2[j-1]) {
                dp[i][j] = dp[i-1][j-1] + 1;
            } else {
                dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
            }
        }
    }
    
    // Backtrack to find actual LCS
    int lcs_len = dp[len1][len2];
    char* lcs = (char*)malloc((lcs_len + 1) * sizeof(char));
    lcs[lcs_len] = '\0';
    
    int i = len1, j = len2;
    int index = lcs_len - 1;
    
    while (i > 0 && j > 0) {
        if (s1[i-1] == s2[j-1]) {
            lcs[index--] = s1[i-1];
            i--;
            j--;
        } else if (dp[i-1][j] > dp[i][j-1]) {
            i--;
        } else {
            j--;
        }
    }
    
    // Free memory
    for (int k = 0; k <= len1; k++) {
        free(dp[k]);
    }
    free(dp);
    
    return lcs;
}

int main() {
    char s1[MAX_LEN], s2[MAX_LEN];
    
    // Read input strings (assuming they are on separate lines)
    printf("Enter first DNA string: ");
    fgets(s1, MAX_LEN, stdin);
    s1[strcspn(s1, "\n")] = 0;  // Remove newline
    
    printf("Enter second DNA string: ");
    fgets(s2, MAX_LEN, stdin);
    s2[strcspn(s2, "\n")] = 0;  // Remove newline
    
    // Find and print LCS
    char* result = findLCS(s1, s2);
    printf("Longest Common Subsequence: %s\n", result);
    
    free(result);
    return 0;
}
```

## Example Usage

**Input:**
```
Enter first DNA string: ATAGGTCG
Enter second DNA string: TATATCGT
```

**Output:**
```
Longest Common Subsequence: TACGT
```

## Time and Space Complexity
- **Time Complexity:** O(m × n) where m and n are lengths of the two strings
- **Space Complexity:** O(m × n) for the DP table

## How it Works
1. The algorithm builds a 2D table where each cell `dp[i][j]` represents the length of LCS for substrings s1[0...i-1] and s2[0...j-1]
2. When characters match, we add 1 to the diagonal value
3. When they don't match, we take the maximum of the left or top cell
4. Backtracking reconstructs the actual subsequence by following the path from dp[m][n] back to dp[0][0]

This solution handles the standard dynamic programming approach for finding the longest common subsequence efficiently.