# Rosalind Problem: Counting Point Mutations (BA1G)

## Problem Statement
Given two strings s and t of equal length, count the number of positions where the characters differ (Hamming distance).

## Solution

```c
#include <stdio.h>
#include <string.h>

int main() {
    char s[1000], t[1000];
    
    // Read the two DNA strings
    scanf("%s", s);
    scanf("%s", t);
    
    int length = strlen(s);
    int mutations = 0;
    
    // Count positions where characters differ
    for (int i = 0; i < length; i++) {
        if (s[i] != t[i]) {
            mutations++;
        }
    }
    
    printf("%d\n", mutations);
    
    return 0;
}
```

## Explanation

1. **Input Reading**: We read two DNA strings `s` and `t` using `scanf()`
2. **Length Calculation**: We determine the length of the strings using `strlen()`
3. **Mutation Counting**: We iterate through each position in both strings and compare characters
4. **Output**: Print the total count of differing positions

## Example

**Input:**
```
GAGCCTACTAACGGGAT
CATCGGTATCGCCAT
```

**Output:**
```
7
```

## Time and Space Complexity

- **Time Complexity**: O(n) where n is the length of the strings
- **Space Complexity**: O(n) for storing the two strings

## Alternative Implementation with File Input

```c
#include <stdio.h>
#include <string.h>

int main() {
    FILE *file = fopen("rosalind_hamm.txt", "r");
    char s[1000], t[1000];
    
    fscanf(file, "%s", s);
    fscanf(file, "%s", t);
    
    int length = strlen(s);
    int mutations = 0;
    
    for (int i = 0; i < length; i++) {
        if (s[i] != t[i]) {
            mutations++;
        }
    }
    
    printf("%d\n", mutations);
    
    fclose(file);
    return 0;
}
```