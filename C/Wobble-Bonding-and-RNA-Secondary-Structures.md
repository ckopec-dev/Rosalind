# Rosalind Problem: Wobble Bonding and RNA Secondary Structures

This problem asks us to count the number of possible secondary structures for an RNA sequence, considering wobble base pairing (G-U pairs are allowed).

## Problem Understanding

We need to count the number of ways to form valid secondary structures where:
- Each nucleotide can pair with another nucleotide
- Valid pairs: A-U, U-A, C-G, G-C, G-U, U-G
- No crossing bonds (non-crossing matching)
- The total count is modulo 10^9 + 7

## Approach

This is a classic dynamic programming problem. We'll use the **Nussinov algorithm** approach with wobble pairing.

Let `dp[i][j]` = number of valid secondary structures for substring from index i to j.

Base cases:
- `dp[i][i] = 1` (single nucleotide can't pair, so 1 way)
- `dp[i][i+1] = 1` (two nucleotides can either pair or not)

Recurrence relation:
- If `s[i]` and `s[j]` can pair: `dp[i][j] += dp[i+1][j-1] + sum(dp[i][k] * dp[k+1][j])`
- Else: `dp[i][j] = sum(dp[i][k] * dp[k+1][j])`

## Solution in C

```c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MOD 1000000007
#define MAXN 505

long long dp[MAXN][MAXN];
char sequence[MAXN];

// Check if two nucleotides can pair (including wobble pairs)
int can_pair(char a, char b) {
    if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A'))
        return 1;
    if ((a == 'C' && b == 'G') || (a == 'G' && b == 'C'))
        return 1;
    if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G'))
        return 1;
    return 0;
}

int main() {
    // Read the RNA sequence
    scanf("%s", sequence);
    int n = strlen(sequence);
    
    // Initialize dp table
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dp[i][j] = 0;
        }
    }
    
    // Base case: single nucleotides
    for (int i = 0; i < n; i++) {
        dp[i][i] = 1;
    }
    
    // Fill the DP table
    for (int len = 2; len <= n; len++) {
        for (int i = 0; i <= n - len; i++) {
            int j = i + len - 1;
            
            // Case 1: No pairing at positions i and j
            dp[i][j] = dp[i+1][j];
            
            // Case 2: Pairing at positions i and j (if possible)
            if (can_pair(sequence[i], sequence[j])) {
                dp[i][j] = (dp[i][j] + dp[i+1][j-1]) % MOD;
            }
            
            // Case 3: Split the structure
            for (int k = i+1; k < j; k++) {
                if (can_pair(sequence[i], sequence[j])) {
                    dp[i][j] = (dp[i][j] + (dp[i+1][k] * dp[k+1][j]) % MOD) % MOD;
                }
            }
        }
    }
    
    printf("%lld\n", dp[0][n-1]);
    
    return 0;
}
```

Wait, I need to correct this approach. The recurrence is more complex for RNA secondary structure prediction.

```c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MOD 1000000007
#define MAXN 505

long long dp[MAXN][MAXN];
char sequence[MAXN];

// Check if two nucleotides can pair (including wobble pairs)
int can_pair(char a, char b) {
    if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A'))
        return 1;
    if ((a == 'C' && b == 'G') || (a == 'G' && b == 'C'))
        return 1;
    if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G'))
        return 1;
    return 0;
}

int main() {
    // Read the RNA sequence
    scanf("%s", sequence);
    int n = strlen(sequence);
    
    // Initialize dp table
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dp[i][j] = 0;
        }
    }
    
    // Base case: single nucleotides
    for (int i = 0; i < n; i++) {
        dp[i][i] = 1;
    }
    
    // Fill the DP table using correct recurrence
    for (int len = 2; len <= n; len++) {
        for (int i = 0; i <= n - len; i++) {
            int j = i + len - 1;
            
            // No pairing at ends
            dp[i][j] = dp[i+1][j];
            
            // If we can pair the ends, add this case
            if (can_pair(sequence[i], sequence[j])) {
                dp[i][j] = (dp[i][j] + dp[i+1][j-1]) % MOD;
            }
            
            // Try all possible splits
            for (int k = i+1; k < j; k++) {
                if (can_pair(sequence[i], sequence[j])) {
                    dp[i][j] = (dp[i][j] + (dp[i][k] * dp[k+1][j]) % MOD) % MOD;
                }
            }
        }
    }
    
    printf("%lld\n", dp[0][n-1]);
    
    return 0;
}
```

Actually, let me re-read the problem and provide a correct solution:

```c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MOD 1000000007
#define MAXN 505

long long dp[MAXN][MAXN];
char sequence[MAXN];

// Check if two nucleotides can pair (including wobble pairs)
int can_pair(char a, char b) {
    if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A'))
        return 1;
    if ((a == 'C' && b == 'G') || (a == 'G' && b == 'C'))
        return 1;
    if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G'))
        return 1;
    return 0;
}

int main() {
    // Read the RNA sequence
    scanf("%s", sequence);
    int n = strlen(sequence);
    
    // Initialize dp table
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dp[i][j] = 0;
        }
    }
    
    // Base case: single nucleotides
    for (int i = 0; i < n; i++) {
        dp[i][i] = 1;
    }
    
    // Fill the DP table using correct recurrence
    for (int len = 2; len <= n; len++) {
        for (int i = 0; i <= n - len; i++) {
            int j = i + len - 1;
            
            // Case 1: No pairing at positions i and j
            dp[i][j] = dp[i+1][j];
            
            // Case 2: Pairing at positions i and j (if possible)
            if (can_pair(sequence[i], sequence[j])) {
                if (i+1 <= j-1) {
                    dp[i][j] = (dp[i][j] + dp[i+1][j-1]) % MOD;
                } else {
                    dp[i][j] = (dp[i][j] + 1) % MOD;
                }
            }
            
            // Case 3: Split the structure
            for (int k = i+1; k < j; k++) {
                dp[i][j] = (dp[i][j] + (dp[i][k] * dp[k+1][j]) % MOD) % MOD;
            }
        }
    }
    
    printf("%lld\n", dp[0][n-1]);
    
    return 0;
}
```

This solution implements the Nussinov algorithm for RNA secondary structure prediction with wobble base pairing. The time complexity is O(n³) and space complexity is O(n²).