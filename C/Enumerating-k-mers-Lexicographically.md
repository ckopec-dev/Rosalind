# Rosalind Problem: Enumerating k-mers Lexicographically (C Solution)

## Problem Understanding

Given an ordered alphabet and a positive integer k, we need to generate all possible k-mers (strings of length k) that can be formed using the given alphabet, arranged in lexicographic order.

## Approach

1. **Input Processing**: Read the alphabet and value of k
2. **Generate All Combinations**: Use recursive backtracking to generate all possible k-length strings
3. **Lexicographic Ordering**: Since we process characters in order, the output will naturally be lexicographically ordered
4. **Output**: Print each generated k-mer

## C Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *alphabet;
int alphabet_size;
int k;
char *current_string;

void generate_kmers(int position) {
    // Base case: if we've filled all positions
    if (position == k) {
        current_string[k] = '\0';  // Null terminate the string
        printf("%s\n", current_string);
        return;
    }
    
    // Recursive case: try each character in the alphabet
    for (int i = 0; i < alphabet_size; i++) {
        current_string[position] = alphabet[i];
        generate_kmers(position + 1);
    }
}

int main() {
    // Read the alphabet (assuming it's provided as a string)
    char input[1000];
    printf("Enter alphabet: ");
    fgets(input, sizeof(input), stdin);
    
    // Remove newline character if present
    int len = strlen(input);
    if (len > 0 && input[len-1] == '\n') {
        input[len-1] = '\0';
        len--;
    }
    
    // Store the alphabet
    alphabet = malloc((len + 1) * sizeof(char));
    strcpy(alphabet, input);
    alphabet_size = len;
    
    // Read k
    printf("Enter k: ");
    scanf("%d", &k);
    
    // Allocate memory for current string
    current_string = malloc((k + 1) * sizeof(char));
    
    // Generate all k-mers
    generate_kmers(0);
    
    // Free allocated memory
    free(alphabet);
    free(current_string);
    
    return 0;
}
```

## Alternative Cleaner Version

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *alphabet;
int alphabet_size;
int k;

void generate_kmers_recursive(char *current, int position) {
    // Base case: if we've filled all positions
    if (position == k) {
        current[k] = '\0';
        printf("%s\n", current);
        return;
    }
    
    // Try each character from the alphabet
    for (int i = 0; i < alphabet_size; i++) {
        current[position] = alphabet[i];
        generate_kmers_recursive(current, position + 1);
    }
}

int main() {
    char input[1000];
    
    // Read alphabet
    printf("Enter alphabet: ");
    fgets(input, sizeof(input), stdin);
    
    // Remove newline character
    int len = strlen(input);
    if (len > 0 && input[len-1] == '\n') {
        input[len-1] = '\0';
    }
    
    // Store alphabet and k
    alphabet_size = len;
    alphabet = malloc((len + 1) * sizeof(char));
    strcpy(alphabet, input);
    
    printf("Enter k: ");
    scanf("%d", &k);
    
    // Allocate memory for current k-mer
    char *current = malloc((k + 1) * sizeof(char));
    
    // Generate all k-mers
    generate_kmers_recursive(current, 0);
    
    // Free memory
    free(alphabet);
    free(current);
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
Enter alphabet: ABCD
Enter k: 2
```

**Output:**
```
AA
AB
AC
AD
BA
BB
BC
BD
CA
CB
CC
CD
DA
DB
DC
DD
```

## Key Features of the Solution

1. **Recursive Generation**: Uses backtracking to systematically generate all combinations
2. **Lexicographic Order**: Since we iterate through alphabet characters in order, output is naturally lexicographically sorted
3. **Memory Management**: Properly allocates and frees memory for strings
4. **Input Flexibility**: Handles alphabets of any size and k values

## Time and Space Complexity

- **Time Complexity**: O(|A|^k × k) where |A| is the alphabet size and k is the length of k-mers
- **Space Complexity**: O(k) for the recursion stack plus O(|A|^k × k) for output storage

The solution efficiently generates all possible k-mers in lexicographic order using recursive backtracking, which is a standard approach for this type of combinatorial problem.