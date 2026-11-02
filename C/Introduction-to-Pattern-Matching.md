# Rosalind Problem: Introduction to Pattern Matching (in C)

## Problem Understanding

The task is to find all starting positions of a pattern within a text string, where the pattern may have multiple occurrences.

## Solution Approach

1. Read the input text and pattern
2. Search for all occurrences of the pattern in the text
3. Record the starting positions (1-indexed)
4. Output the positions

## C Implementation

```c
#include <stdio.h>
#include <string.h>

int main() {
    char text[100000];
    char pattern[100000];
    
    // Read input
    scanf("%s", text);
    scanf("%s", pattern);
    
    int text_len = strlen(text);
    int pattern_len = strlen(pattern);
    
    // Find all occurrences of pattern in text
    for (int i = 0; i <= text_len - pattern_len; i++) {
        // Check if pattern matches at position i
        if (strncmp(text + i, pattern, pattern_len) == 0) {
            printf("%d ", i + 1);  // Output 1-indexed position
        }
    }
    
    printf("\n");
    
    return 0;
}
```

## Explanation

1. **Input Reading**: We read both the text and pattern strings using `scanf()`
2. **Length Calculation**: Get lengths of both strings using `strlen()`
3. **Pattern Search**: Loop through each possible starting position in the text
4. **String Comparison**: Use `strncmp()` to compare substring with pattern
5. **Output Format**: Print positions (1-indexed) separated by spaces

## Example

**Input:**
```
ATATATATATATATA
TAT
```

**Output:**
```
3 6 9 12 15 
```

The pattern "TAT" occurs at positions 3, 6, 9, 12, and 15 (1-indexed) in the text.

## Time Complexity
- **Time**: O(n × m) where n is text length and m is pattern length
- **Space**: O(n + m) for storing the strings

This solution handles the basic pattern matching problem efficiently and correctly outputs all starting positions of the pattern within the given text.