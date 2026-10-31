# Rosalind Problem: Edit Distance (Edit_Distance)

## Problem Description
The edit distance between two strings is the minimum number of operations (insertions, deletions, or substitutions) required to transform one string into another.

## Solution Approach
We'll use dynamic programming to solve this problem. The DP table `dp[i][j]` represents the minimum edit distance between the first `i` characters of string1 and the first `j` characters of string2.

## C Implementation

```c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int min(int a, int b, int c) {
    int min_val = a;
    if (b < min_val) min_val = b;
    if (c < min_val) min_val = c;
    return min_val;
}

int edit_distance(char *str1, char *str2) {
    int m = strlen(str1);
    int n = strlen(str2);
    
    // Create DP table
    int **dp = (int**)malloc((m + 1) * sizeof(int*));
    for (int i = 0; i <= m; i++) {
        dp[i] = (int*)malloc((n + 1) * sizeof(int));
    }
    
    // Initialize base cases
    // Empty string to str2[0:j] requires j insertions
    for (int j = 0; j <= n; j++) {
        dp[0][j] = j;
    }
    
    // str1[0:i] to empty string requires i deletions
    for (int i = 0; i <= m; i++) {
        dp[i][0] = i;
    }
    
    // Fill the DP table
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (str1[i-1] == str2[j-1]) {
                // Characters match, no operation needed
                dp[i][j] = dp[i-1][j-1];
            } else {
                // Take minimum of three operations
                dp[i][j] = 1 + min(dp[i-1][j],     // deletion
                                   dp[i][j-1],     // insertion
                                   dp[i-1][j-1]);  // substitution
            }
        }
    }
    
    int result = dp[m][n];
    
    // Free memory
    for (int i = 0; i <= m; i++) {
        free(dp[i]);
    }
    free(dp);
    
    return result;
}

int main() {
    char str1[1000], str2[1000];
    
    printf("Enter first string: ");
    scanf("%s", str1);
    printf("Enter second string: ");
    scanf("%s", str2);
    
    int distance = edit_distance(str1, str2);
    printf("Edit distance: %d\n", distance);
    
    return 0;
}
```

## Example Usage

**Input:**
```
Enter first string: AGTCTGAT
Enter second string: AGTAATGT
```

**Output:**
```
Edit distance: 2
```

## Algorithm Explanation

1. **DP Table Construction**: Create a table `dp[m+1][n+1]` where `m` and `n` are the lengths of the two strings.

2. **Base Cases**: 
   - `dp[0][j] = j`: Converting empty string to first `j` characters requires `j` insertions
   - `dp[i][0] = i`: Converting first `i` characters to empty string requires `i` deletions

3. **Recurrence Relation**:
   - If `str1[i-1] == str2[j-1]`: No operation needed, `dp[i][j] = dp[i-1][j-1]`
   - Otherwise: Take minimum of three operations plus 1, `dp[i][j] = 1 + min(dp[i-1][j], dp[i][j-1], dp[i-1][j-1])`

4. **Result**: The value at `dp[m][n]` gives the edit distance.

## Time and Space Complexity
- **Time Complexity**: O(m × n) where m and n are the lengths of the strings
- **Space Complexity**: O(m × n) for the DP table

This solution efficiently computes the minimum number of operations needed to transform one string into another using dynamic programming.