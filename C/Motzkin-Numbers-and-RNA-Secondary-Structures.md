# Rosalind Problem: Motzkin Numbers and RNA Secondary Structures

## Problem Understanding

This problem asks us to compute the number of possible secondary structures for an RNA sequence, which is related to Motzkin numbers. The key constraints are:
- Each nucleotide can pair with its complement (A with U, G with C)
- Pairs cannot cross each other (no pseudoknots)
- Each nucleotide can participate in at most one pair

This is equivalent to counting the number of non-crossing partitions of a set, which relates to Catalan numbers and Motzkin numbers.

## Approach

We'll use dynamic programming where `dp[i][j]` represents the number of valid secondary structures for substring from index `i` to `j`.

The recurrence relation:
1. If we don't pair positions `i` and `j`: `dp[i][j] += dp[i+1][j]`
2. If we pair position `i` with some position `k` (where `k` can be any valid pairing position): 
   `dp[i][j] += dp[i+1][k-1] * dp[k+1][j]`

## Solution

```c
#include <stdio.h>
#include <string.h>

#define MAX_N 300
#define MOD 1000000

long long dp[MAX_N][MAX_N];

// Check if two nucleotides can pair
int can_pair(char a, char b) {
    return (a == 'A' && b == 'U') || 
           (a == 'U' && b == 'A') ||
           (a == 'G' && b == 'C') ||
           (a == 'C' && b == 'G');
}

int main() {
    char sequence[MAX_N];
    
    // Read the RNA sequence
    scanf("%s", sequence);
    int n = strlen(sequence);
    
    // Initialize dp table
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dp[i][j] = 0;
        }
    }
    
    // Base cases: empty substring has one structure (empty structure)
    for (int i = 0; i < n; i++) {
        dp[i][i] = 1;
    }
    
    // Fill the dp table
    for (int len = 2; len <= n; len++) {
        for (int i = 0; i <= n - len; i++) {
            int j = i + len - 1;
            
            // Case 1: don't pair the first nucleotide
            dp[i][j] += dp[i+1][j];
            dp[i][j] %= MOD;
            
            // Case 2: pair the first nucleotide with some k
            for (int k = i + 1; k <= j; k++) {
                if (can_pair(sequence[i], sequence[k])) {
                    long long temp = (dp[i+1][k-1] * dp[k+1][j]) % MOD;
                    dp[i][j] += temp;
                    dp[i][j] %= MOD;
                }
            }
        }
    }
    
    printf("%lld\n", dp[0][n-1]);
    
    return 0;
}
```

## Explanation

1. **Input Processing**: Read the RNA sequence from standard input.

2. **Dynamic Programming Setup**: 
   - `dp[i][j]` represents the number of valid secondary structures for substring from index `i` to `j`
   - Initialize base cases where single nucleotides have one structure (empty structure)

3. **Recurrence Relation**:
   - For each substring, we consider two cases:
     1. Don't pair the first nucleotide: `dp[i+1][j]`
     2. Pair the first nucleotide with a valid partner `k`: `dp[i+1][k-1] * dp[k+1][j]`
   - A valid pairing means the nucleotides are complementary (A-U, U-A, G-C, C-G)

4. **Pairing Validation**: The `can_pair` function checks if two nucleotides can form a base pair.

5. **Modular Arithmetic**: Results are taken modulo 1000000 to prevent overflow.

## Time and Space Complexity

- **Time Complexity**: O(n³) where n is the length of the RNA sequence
- **Space Complexity**: O(n²) for the DP table

This solution handles the constraint that RNA secondary structures must be non-crossing, which makes it equivalent to counting valid bracketings or non-crossing partitions.