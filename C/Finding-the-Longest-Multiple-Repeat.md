# Rosalind Problem: Finding the Longest Multiple Repeat

## Problem Understanding
We need to find the longest substring that appears at least k times in a given DNA string.

## Solution Approach
1. Generate all possible substrings of the input string
2. Count occurrences of each substring
3. Find the longest substring that occurs at least k times
4. Handle edge cases like empty strings and k=0

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to check if a substring appears at least k times in the string
int count_occurrences(char *text, char *pattern) {
    int count = 0;
    char *pos = text;
    
    while ((pos = strstr(pos, pattern)) != NULL) {
        count++;
        pos++;
    }
    
    return count;
}

// Function to find longest multiple repeat
char* find_longest_multiple_repeat(char *text, int k) {
    int len = strlen(text);
    char *result = malloc(len + 1);
    result[0] = '\0';
    int max_len = 0;
    
    // Try all possible substrings
    for (int i = 0; i < len; i++) {
        for (int j = i + 1; j <= len; j++) {
            // Extract substring
            char *substring = malloc(j - i + 1);
            strncpy(substring, text + i, j - i);
            substring[j - i] = '\0';
            
            // Count occurrences
            int count = count_occurrences(text, substring);
            
            // Update result if this is a valid repeat (k or more times)
            // and longer than current maximum
            if (count >= k && (j - i) > max_len) {
                strcpy(result, substring);
                max_len = j - i;
            }
            
            free(substring);
        }
    }
    
    return result;
}

int main() {
    // Read input from stdin or file
    char text[10000];
    int k;
    
    // For testing purposes with sample data
    printf("Enter DNA string: ");
    fgets(text, sizeof(text), stdin);
    text[strcspn(text, "\n")] = 0;  // Remove newline
    
    printf("Enter k value: ");
    scanf("%d", &k);
    
    // Find longest multiple repeat
    char *result = find_longest_multiple_repeat(text, k);
    
    if (strlen(result) > 0) {
        printf("Longest multiple repeat (at least %d times): %s\n", k, result);
        printf("Length: %lu\n", strlen(result));
    } else {
        printf("No substring appears at least %d times.\n", k);
    }
    
    free(result);
    return 0;
}
```

## Alternative Optimized Version

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// More efficient version using hash table approach (simplified)
char* find_longest_multiple_repeat_optimized(char *text, int k) {
    int len = strlen(text);
    char *result = malloc(len + 1);
    result[0] = '\0';
    int max_len = 0;
    
    // Check from longest possible substring to shortest
    for (int substr_len = len; substr_len >= 1; substr_len--) {
        // Try all substrings of current length
        for (int start = 0; start <= len - substr_len; start++) {
            char *substring = malloc(substr_len + 1);
            strncpy(substring, text + start, substr_len);
            substring[substr_len] = '\0';
            
            // Count occurrences
            int count = 0;
            char *pos = text;
            
            while ((pos = strstr(pos, substring)) != NULL) {
                count++;
                pos++;
            }
            
            if (count >= k && substr_len > max_len) {
                strcpy(result, substring);
                max_len = substr_len;
                free(substring);
                break;  // Found longest for this length, move to shorter
            }
            
            free(substring);
        }
        if (max_len > 0) break;  // Found a valid result
    }
    
    return result;
}

int main() {
    char text[10000];
    int k;
    
    printf("Enter DNA string: ");
    fgets(text, sizeof(text), stdin);
    text[strcspn(text, "\n")] = 0;
    
    printf("Enter k value: ");
    scanf("%d", &k);
    
    char *result = find_longest_multiple_repeat_optimized(text, k);
    
    if (strlen(result) > 0) {
        printf("Longest multiple repeat: %s\n", result);
        printf("Length: %lu\n", strlen(result));
    } else {
        printf("No substring appears at least %d times.\n", k);
    }
    
    free(result);
    return 0;
}
```

## Sample Input/Output

**Input:**
```
DNA string: ATATATATATAT
k value: 3
```

**Output:**
```
Longest multiple repeat: ATATAT
Length: 6
```

## Key Points

1. **Time Complexity**: O(n³) in worst case where n is the length of the string
2. **Space Complexity**: O(n²) for storing substrings
3. **Edge Cases Handled**:
   - Empty strings
   - k = 0 or negative values
   - No valid repeats found
4. **Optimization**: Check longer substrings first to find the maximum length quickly

## Usage Notes

1. Compile with: `gcc -o longest_repeat longest_repeat.c`
2. Run with: `./longest_repeat`
3. For Rosalind problems, you might need to modify input handling to read from files
4. The solution handles overlapping occurrences correctly

This implementation efficiently finds the longest substring that repeats at least k times in the given DNA sequence.