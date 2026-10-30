# Rosalind Problem: Ordering Strings of Varying Length Lexicographically

## Problem Statement
Given a collection of strings, return all permutations of the strings in lexicographic order.

## Solution Approach
1. Read input strings from stdin
2. Sort strings lexicographically using qsort
3. Print sorted strings

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Comparison function for qsort
int compare_strings(const void *a, const void *b) {
    const char **str1 = (const char **)a;
    const char **str2 = (const char **)b;
    return strcmp(*str1, *str2);
}

int main() {
    // Read number of strings
    int n;
    scanf("%d", &n);
    
    // Allocate memory for array of string pointers
    char **strings = malloc(n * sizeof(char *));
    
    // Read each string
    for (int i = 0; i < n; i++) {
        // Assuming maximum length of string is 100
        strings[i] = malloc(101 * sizeof(char));
        scanf("%s", strings[i]);
    }
    
    // Sort strings lexicographically
    qsort(strings, n, sizeof(char *), compare_strings);
    
    // Print sorted strings
    for (int i = 0; i < n; i++) {
        printf("%s\n", strings[i]);
        free(strings[i]);  // Free individual string
    }
    
    free(strings);  // Free array of pointers
    
    return 0;
}
```

## Alternative Solution with Dynamic Input Reading

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Comparison function for qsort
int compare_strings(const void *a, const void *b) {
    const char **str1 = (const char **)a;
    const char **str2 = (const char **)b;
    return strcmp(*str1, *str2);
}

int main() {
    // Read strings until EOF
    char **strings = NULL;
    int capacity = 10;
    int count = 0;
    
    // Allocate initial memory
    strings = malloc(capacity * sizeof(char *));
    
    // Read input strings
    char buffer[101];  // Buffer for reading strings
    while (scanf("%s", buffer) == 1) {
        // Resize array if needed
        if (count >= capacity) {
            capacity *= 2;
            strings = realloc(strings, capacity * sizeof(char *));
        }
        
        // Allocate memory and copy string
        strings[count] = malloc(strlen(buffer) + 1);
        strcpy(strings[count], buffer);
        count++;
    }
    
    // Sort strings lexicographically
    qsort(strings, count, sizeof(char *), compare_strings);
    
    // Print sorted strings
    for (int i = 0; i < count; i++) {
        printf("%s\n", strings[i]);
        free(strings[i]);  // Free individual string
    }
    
    free(strings);  // Free array of pointers
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
6
AAC
AAT
ACT
AGT
CAT
GAT
```

**Output:**
```
AAC
AAT
ACT
AGT
CAT
GAT
```

## Key Points

1. **Memory Management**: Properly allocate and free memory for strings
2. **String Comparison**: Use `strcmp()` for lexicographic comparison
3. **Sorting**: Use `qsort()` with custom comparison function
4. **Input Handling**: Handle variable number of input strings
5. **Error Checking**: The second solution handles dynamic input size

## Time Complexity
- O(n log n) for sorting where n is the number of strings
- O(m × n) for string comparisons where m is average string length

## Space Complexity
- O(n × m) for storing all strings