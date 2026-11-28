# Rosalind Problem: Compute the Size of a Spectral Dictionary

## Problem Understanding

The problem asks us to compute the size of a spectral dictionary, which is essentially counting how many distinct k-mers (substrings of length k) can be formed from a given string.

## Approach

1. Read the input string and k value
2. Generate all possible substrings of length k
3. Count unique substrings using a set-like data structure
4. Return the count

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Hash function for strings
unsigned int hash(const char* str, int len) {
    unsigned int hash = 5381;
    for (int i = 0; i < len; i++) {
        hash = ((hash << 5) + hash) + str[i];
    }
    return hash;
}

// Simple hash set implementation
typedef struct Node {
    char* key;
    struct Node* next;
} Node;

typedef struct {
    Node** buckets;
    int size;
    int capacity;
} HashSet;

// Initialize hash set
HashSet* create_hashset(int capacity) {
    HashSet* set = (HashSet*)malloc(sizeof(HashSet));
    set->capacity = capacity;
    set->size = 0;
    set->buckets = (Node**)calloc(capacity, sizeof(Node*));
    return set;
}

// Hash function for strings
int hash_string(const char* str, int len, int capacity) {
    unsigned int hash_value = 5381;
    for (int i = 0; i < len; i++) {
        hash_value = ((hash_value << 5) + hash_value) + str[i];
    }
    return hash_value % capacity;
}

// Check if string exists in hash set
int contains(HashSet* set, const char* str, int len) {
    int index = hash_string(str, len, set->capacity);
    Node* current = set->buckets[index];
    
    while (current != NULL) {
        if (strlen(current->key) == len && strncmp(current->key, str, len) == 0) {
            return 1;
        }
        current = current->next;
    }
    return 0;
}

// Insert string into hash set
void insert(HashSet* set, const char* str, int len) {
    if (contains(set, str, len)) {
        return; // Already exists
    }
    
    int index = hash_string(str, len, set->capacity);
    Node* new_node = (Node*)malloc(sizeof(Node));
    new_node->key = (char*)malloc(len + 1);
    strncpy(new_node->key, str, len);
    new_node->key[len] = '\0';
    new_node->next = set->buckets[index];
    set->buckets[index] = new_node;
    set->size++;
}

int main() {
    char text[10000];
    int k;
    
    // Read input
    scanf("%s", text);
    scanf("%d", &k);
    
    int n = strlen(text);
    HashSet* dict = create_hashset(10000);
    
    // Generate all substrings of length k and add to set
    for (int i = 0; i <= n - k; i++) {
        char substring[1000];
        strncpy(substring, text + i, k);
        substring[k] = '\0';
        insert(dict, substring, k);
    }
    
    printf("%d\n", dict->size);
    
    // Free memory
    for (int i = 0; i < dict->capacity; i++) {
        Node* current = dict->buckets[i];
        while (current != NULL) {
            Node* temp = current;
            current = current->next;
            free(temp->key);
            free(temp);
        }
    }
    free(dict->buckets);
    free(dict);
    
    return 0;
}
```

## Explanation

1. **Input Reading**: The program reads the input string and k value from standard input.

2. **Hash Set Implementation**: 
   - Uses a simple hash table with chaining to store unique substrings
   - Hash function uses polynomial rolling hash technique
   - Handles collisions through linked lists

3. **Substring Generation**:
   - Iterates through all possible starting positions (0 to n-k)
   - Extracts substrings of length k
   - Inserts each substring into the hash set

4. **Counting**: 
   - The hash set automatically handles duplicates
   - Returns the final size of the set

## Time Complexity
- O(n×k) where n is the length of the string and k is the substring length
- In practice, hash collisions make it faster than O(n×k) for most cases

## Space Complexity
- O(n×k) in worst case for storing all unique substrings

This solution efficiently handles the spectral dictionary problem by counting unique k-length substrings in a given string.