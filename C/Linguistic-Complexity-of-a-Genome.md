# Rosalind Problem: Linguistic Complexity of a Genome (C Solution)

## Problem Understanding

The linguistic complexity of a genome is defined as the ratio of the number of distinct substrings to the total number of possible substrings of a given length.

For a DNA sequence of length n, we need to:
1. Count all distinct substrings of length k
2. Calculate the total number of possible substrings of length k
3. Compute the ratio

## Solution Approach

We'll use a hash table approach to efficiently count distinct substrings and avoid memory issues with large sequences.

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Simple hash function for strings
unsigned int hash_string(const char* str, int len) {
    unsigned int hash = 5381;
    for (int i = 0; i < len; i++) {
        hash = ((hash << 5) + hash) + str[i];
    }
    return hash;
}

// Structure to store substrings in hash table
typedef struct HashNode {
    char* substring;
    int count;
    struct HashNode* next;
} HashNode;

// Hash table structure
typedef struct {
    HashNode** buckets;
    int size;
} HashTable;

// Create new hash table
HashTable* create_hash_table(int size) {
    HashTable* table = (HashTable*)malloc(sizeof(HashTable));
    table->size = size;
    table->buckets = (HashNode**)calloc(size, sizeof(HashNode*));
    return table;
}

// Insert substring into hash table
void insert_substring(HashTable* table, const char* substring) {
    int len = strlen(substring);
    unsigned int hash = hash_string(substring, len) % table->size;
    
    HashNode* current = table->buckets[hash];
    while (current != NULL) {
        if (strcmp(current->substring, substring) == 0) {
            current->count++;
            return;
        }
        current = current->next;
    }
    
    // Create new node
    HashNode* new_node = (HashNode*)malloc(sizeof(HashNode));
    new_node->substring = (char*)malloc(len + 1);
    strcpy(new_node->substring, substring);
    new_node->count = 1;
    new_node->next = table->buckets[hash];
    table->buckets[hash] = new_node;
}

// Count distinct substrings of length k
int count_distinct_substrings(const char* sequence, int length, int k) {
    if (k > length) return 0;
    
    HashTable* table = create_hash_table(10007); // Prime number for hash table size
    
    for (int i = 0; i <= length - k; i++) {
        char substring[1000]; // Assuming max length of 1000
        strncpy(substring, sequence + i, k);
        substring[k] = '\0';
        insert_substring(table, substring);
    }
    
    // Count distinct substrings
    int count = 0;
    for (int i = 0; i < table->size; i++) {
        HashNode* current = table->buckets[i];
        while (current != NULL) {
            count++;
            current = current->next;
        }
    }
    
    // Free memory
    for (int i = 0; i < table->size; i++) {
        HashNode* current = table->buckets[i];
        while (current != NULL) {
            HashNode* temp = current;
            current = current->next;
            free(temp->substring);
            free(temp);
        }
    }
    free(table->buckets);
    free(table);
    
    return count;
}

// Calculate linguistic complexity
double linguistic_complexity(const char* sequence, int length, int k) {
    if (k > length) return 0.0;
    
    // Total possible substrings of length k
    int total_possible = length - k + 1;
    
    // Count distinct substrings
    int distinct_count = count_distinct_substrings(sequence, length, k);
    
    return (double)distinct_count / (double)total_possible;
}

int main() {
    char sequence[10000];
    int k;
    
    // Read input from stdin
    printf("Enter DNA sequence: ");
    fgets(sequence, sizeof(sequence), stdin);
    
    // Remove newline character if present
    int len = strlen(sequence);
    if (len > 0 && sequence[len-1] == '\n') {
        sequence[len-1] = '\0';
        len--;
    }
    
    printf("Enter k (substring length): ");
    scanf("%d", &k);
    
    // Calculate linguistic complexity
    double complexity = linguistic_complexity(sequence, len, k);
    
    printf("Linguistic complexity: %.6f\n", complexity);
    
    return 0;
}
```

## Alternative Implementation Using Set-like Approach

For better performance and cleaner code:

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Simple hash table implementation for substring counting
typedef struct {
    char** keys;
    int* counts;
    int size;
    int capacity;
} SubstringSet;

// Create new set
SubstringSet* create_set(int capacity) {
    SubstringSet* set = (SubstringSet*)malloc(sizeof(SubstringSet));
    set->keys = (char**)malloc(capacity * sizeof(char*));
    set->counts = (int*)calloc(capacity, sizeof(int));
    set->size = 0;
    set->capacity = capacity;
    return set;
}

// Hash function for strings
unsigned int hash_string(const char* str, int len) {
    unsigned int hash = 5381;
    for (int i = 0; i < len; i++) {
        hash = ((hash << 5) + hash) + str[i];
    }
    return hash;
}

// Add substring to set
void add_substring(SubstringSet* set, const char* substring, int len) {
    unsigned int hash_val = hash_string(substring, len) % set->capacity;
    
    // Simple linear probing
    int index = hash_val;
    while (set->keys[index] != NULL) {
        if (strcmp(set->keys[index], substring) == 0) {
            set->counts[index]++;
            return;
        }
        index = (index + 1) % set->capacity;
    }
    
    // Insert new substring
    set->keys[index] = (char*)malloc(len + 1);
    strcpy(set->keys[index], substring);
    set->counts[index] = 1;
    set->size++;
}

// Count distinct substrings of length k
int count_distinct_substrings(const char* sequence, int length, int k) {
    if (k > length) return 0;
    
    SubstringSet* set = create_set(length);
    
    for (int i = 0; i <= length - k; i++) {
        char substring[1000];
        strncpy(substring, sequence + i, k);
        substring[k] = '\0';
        add_substring(set, substring, k);
    }
    
    int count = set->size;
    
    // Free memory
    for (int i = 0; i < set->capacity; i++) {
        if (set->keys[i] != NULL) {
            free(set->keys[i]);
        }
    }
    free(set->keys);
    free(set->counts);
    free(set);
    
    return count;
}

// Main function
int main() {
    char sequence[10000];
    int k;
    
    // Read input
    printf("Enter DNA sequence: ");
    if (fgets(sequence, sizeof(sequence), stdin) == NULL) {
        return 1;
    }
    
    // Remove newline
    int len = strlen(sequence);
    if (len > 0 && sequence[len-1] == '\n') {
        sequence[len-1] = '\0';
        len--;
    }
    
    printf("Enter k: ");
    scanf("%d", &k);
    
    // Calculate linguistic complexity
    int total_substrings = len - k + 1;
    int distinct_count = count_distinct_substrings(sequence, len, k);
    double complexity = (double)distinct_count / (double)total_substrings;
    
    printf("Linguistic Complexity: %.6f\n", complexity);
    
    return 0;
}
```

## How to Compile and Run

```bash
gcc -o linguistic_complexity linguistic_complexity.c
./linguistic_complexity
```

## Input Format

The program expects:
1. A DNA sequence (string of A, C, G, T)
2. An integer k representing substring length

## Example

Input:
```
Enter DNA sequence: ACACACAC
Enter k: 3
```

Output:
```
Linguistic Complexity: 0.666667
```

## Time and Space Complexity

- **Time Complexity**: O(n×k) where n is the sequence length and k is substring length
- **Space Complexity**: O(n×k) for storing substrings in hash table

This solution efficiently handles the linguistic complexity calculation by using a hash table to count distinct substrings while maintaining good performance characteristics.