# Rosalind Problem: Counting Optimal Alignments

## Problem Understanding

This problem asks us to count the number of optimal global alignments between two sequences using dynamic programming with traceback counting.

## Solution Approach

I'll implement a solution that:
1. Computes the edit distance using dynamic programming
2. Counts the number of optimal alignments by tracking backtracking paths
3. Uses modular arithmetic to handle large numbers

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MOD 1000000
#define MAXLEN 1001

long long dp[MAXLEN][MAXLEN];
long long count[MAXLEN][MAXLEN];

int main() {
    char seq1[MAXLEN], seq2[MAXLEN];
    
    // Read input sequences
    scanf("%s", seq1);
    scanf("%s", seq2);
    
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    
    // Initialize DP table and count table
    for (int i = 0; i <= len1; i++) {
        dp[i][0] = i;
        count[i][0] = (i == 0) ? 1 : 0;
    }
    
    for (int j = 0; j <= len2; j++) {
        dp[0][j] = j;
        count[0][j] = (j == 0) ? 1 : 0;
    }
    
    // Fill the DP table
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            long long match_cost = (seq1[i-1] == seq2[j-1]) ? 0 : 1;
            
            // Calculate minimum cost
            long long min_cost = dp[i-1][j-1] + match_cost;
            min_cost = (min_cost < dp[i-1][j] + 1) ? min_cost : dp[i-1][j] + 1;
            min_cost = (min_cost < dp[i][j-1] + 1) ? min_cost : dp[i][j-1] + 1;
            
            dp[i][j] = min_cost;
            
            // Count paths
            count[i][j] = 0;
            
            if (dp[i-1][j-1] + match_cost == min_cost) {
                count[i][j] = (count[i][j] + count[i-1][j-1]) % MOD;
            }
            
            if (dp[i-1][j] + 1 == min_cost) {
                count[i][j] = (count[i][j] + count[i-1][j]) % MOD;
            }
            
            if (dp[i][j-1] + 1 == min_cost) {
                count[i][j] = (count[i][j] + count[i][j-1]) % MOD;
            }
        }
    }
    
    // Print the number of optimal alignments
    printf("%lld\n", count[len1][len2]);
    
    return 0;
}
```

## Explanation

### Key Components:

1. **DP Table (`dp[i][j]`)**: Stores the minimum edit distance between prefixes of the two sequences
2. **Count Table (`count[i][j]`)**: Stores the number of ways to achieve the optimal edit distance up to position (i,j)
3. **Base Cases**: Initialize first row and column with insertion/deletion costs
4. **Recurrence Relation**: For each cell, compute minimum cost considering:
   - Match/mismatch (diagonal move)
   - Insertion (up move) 
   - Deletion (left move)

### Algorithm Steps:

1. Read two DNA sequences
2. Initialize base cases for empty string alignments
3. Fill DP table by considering all three possible moves:
   - Diagonal (match/mismatch): `dp[i-1][j-1] + cost`
   - Up (insertion): `dp[i-1][j] + 1`  
   - Left (deletion): `dp[i][j-1] + 1`
4. For counting, accumulate paths from all predecessors that give optimal cost
5. Return the count at the bottom-right corner

### Time Complexity:
- **Time**: O(n×m) where n and m are sequence lengths
- **Space**: O(n×m) for DP and count tables

The solution handles the modular arithmetic to prevent integer overflow, which is crucial when counting paths in large alignments.