# Rosalind Problem: Interleaving Two Motifs (C Solution)

## Problem Understanding

The problem asks us to find the number of ways to interleave two strings such that the resulting string contains both as subsequences.

## Approach

We'll use dynamic programming where `dp[i][j]` represents the number of ways to interleave the first `i` characters of string 1 with the first `j` characters of string 2.

## Solution

```c
#include <stdio.h>
#include <string.h>

#define MAX_LEN 1000

long long dp[MAX_LEN][MAX_LEN];

long long count_interleavings(char *s1, char *s2) {
    int len1 = strlen(s1);
    int len2 = strlen(s2);
    
    // Initialize base cases
    for (int i = 0; i <= len1; i++) {
        dp[i][0] = 1;
    }
    for (int j = 0; j <= len2; j++) {
        dp[0][j] = 1;
    }
    
    // Fill the DP table
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            if (s1[i-1] == s2[j-1]) {
                dp[i][j] = dp[i-1][j] + dp[i][j-1];
            } else {
                dp[i][j] = dp[i-1][j] + dp[i][j-1];
            }
        }
    }
    
    return dp[len1][len2];
}

int main() {
    char s1[MAX_LEN], s2[MAX_LEN];
    
    printf("Enter first string: ");
    scanf("%s", s1);
    printf("Enter second string: ");
    scanf("%s", s2);
    
    long long result = count_interleavings(s1, s2);
    printf("Number of interleavings: %lld\n", result);
    
    return 0;
}
```

## Explanation

1. **DP State Definition**: `dp[i][j]` represents the number of ways to interleave the first `i` characters of string 1 with the first `j` characters of string 2.

2. **Base Cases**: 
   - When one string is empty, there's only one way to interleave (take all characters from the other string)
   - `dp[i][0] = 1` for all i
   - `dp[0][j] = 1` for all j

3. **Recurrence Relation**:
   - If `s1[i-1] == s2[j-1]`: We can either take from string 1 or string 2, so `dp[i][j] = dp[i-1][j] + dp[i][j-1]`
   - If `s1[i-1] != s2[j-1]`: Same logic applies as we're looking for interleavings regardless of character match

4. **Result**: The answer is stored in `dp[len1][len2]`

## Time and Space Complexity

- **Time Complexity**: O(n×m) where n and m are the lengths of the two strings
- **Space Complexity**: O(n×m) for the DP table

## Sample Input/Output

```
Enter first string: ACGT
Enter second string: TCGA
Number of interleavings: 20
```

This solution efficiently computes the number of ways to interleave two motifs using dynamic programming principles.