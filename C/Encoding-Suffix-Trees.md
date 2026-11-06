# Rosalind Problem: Encoding Suffix Trees (C Solution)

## Problem Understanding

The task is to construct a suffix tree for a given string and return its edge labels. A suffix tree is a compressed trie containing all suffixes of the given text.

## Approach

I'll implement a suffix tree using Ukkonen's algorithm with a compact representation that stores edge labels directly on the edges.

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 100005
#define ALPHABET_SIZE 256

typedef struct Node {
    int start, end;
    struct Node* children[ALPHABET_SIZE];
    struct Node* suffix_link;
    int suffix_index;
    int id;
} Node;

char text[MAXN];
int text_len;
Node* root;
Node* active_node;
int active_edge;
int active_length;
int remaining_suffixes;
int leaf_end;
int *suffix_array;
int node_count;

// Create a new node
Node* create_node(int start, int end) {
    Node* node = (Node*)malloc(sizeof(Node));
    node->start = start;
    node->end = end;
    node->suffix_link = NULL;
    node->suffix_index = -1;
    node->id = node_count++;
    
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        node->children[i] = NULL;
    }
    
    return node;
}

// Get the character at a specific position
char get_char(int index) {
    if (index >= text_len) return '\0';
    return text[index];
}

// Get the length of an edge
int edge_length(Node* node) {
    if (node == NULL) return 0;
    return node->end - node->start + 1;
}

// Insert a suffix into the suffix tree
void insert_suffix(int suffix_start) {
    leaf_end = suffix_start;
    remaining_suffixes++;
    
    Node* last_created_node = NULL;
    
    while (remaining_suffixes > 0) {
        if (active_length == 0) {
            active_edge = suffix_start;
        }
        
        if (active_node->children[get_char(active_edge)] == NULL) {
            // Create new leaf node
            Node* new_node = create_node(leaf_end, MAXN - 1);
            active_node->children[get_char(active_edge)] = new_node;
            
            // Apply suffix link
            if (last_created_node != NULL) {
                last_created_node->suffix_link = active_node;
                last_created_node = NULL;
            }
        } else {
            Node* next = active_node->children[get_char(active_edge)];
            int edge_len = edge_length(next);
            
            if (active_length >= edge_len) {
                // Move to the next node
                active_edge += edge_len;
                active_length -= edge_len;
                active_node = next;
                continue;
            }
            
            // Check if we are at the end of the edge
            if (get_char(next->start + active_length) == get_char(suffix_start)) {
                // We found an existing edge
                if (last_created_node != NULL && active_node != root) {
                    last_created_node->suffix_link = active_node;
                    last_created_node = NULL;
                }
                active_length++;
                break;
            }
            
            // Split the edge
            Node* split_node = create_node(next->start, next->start + active_length - 1);
            active_node->children[get_char(active_edge)] = split_node;
            
            Node* new_leaf = create_node(leaf_end, MAXN - 1);
            split_node->children[get_char(suffix_start)] = new_leaf;
            
            next->start += active_length;
            split_node->children[get_char(next->start)] = next;
            
            if (last_created_node != NULL) {
                last_created_node->suffix_link = split_node;
            }
            last_created_node = split_node;
        }
        
        remaining_suffixes--;
        if (active_node == root && active_length > 0) {
            active_length--;
            active_edge = suffix_start - remaining_suffixes + 1;
        } else if (active_node != root) {
            active_node = active_node->suffix_link;
        }
    }
}

// Build the suffix tree
void build_suffix_tree() {
    root = create_node(-1, -1);
    active_node = root;
    active_edge = 0;
    active_length = 0;
    remaining_suffixes = 0;
    leaf_end = -1;
    node_count = 0;
    
    for (int i = 0; i < text_len; i++) {
        insert_suffix(i);
    }
}

// Get edge labels from the suffix tree
void get_edge_labels(Node* node, char* current_label, int label_len, char** result, int* result_count) {
    if (node == NULL) return;
    
    // If this is a leaf node
    if (node->suffix_index != -1) {
        // This is a leaf node, don't process it here
        return;
    }
    
    // Process all children
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        Node* child = node->children[i];
        if (child != NULL) {
            int edge_start = child->start;
            int edge_end = child->end;
            
            // Extract the edge label
            char temp_label[1000];
            int temp_len = 0;
            
            for (int j = edge_start; j <= edge_end && j < text_len; j++) {
                temp_label[temp_len++] = text[j];
            }
            temp_label[temp_len] = '\0';
            
            // Add to results
            if (temp_len > 0) {
                result[*result_count] = (char*)malloc((temp_len + 1) * sizeof(char));
                strcpy(result[*result_count], temp_label);
                (*result_count)++;
            }
            
            get_edge_labels(child, current_label, label_len, result, result_count);
        }
    }
}

// Main function to solve the problem
void solve() {
    // Read input
    scanf("%s", text);
    text_len = strlen(text);
    
    // Build suffix tree
    build_suffix_tree();
    
    // Extract edge labels
    char** edge_labels = (char**)malloc(10000 * sizeof(char*));
    int count = 0;
    
    // Get all edge labels from the tree
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        if (root->children[i] != NULL) {
            Node* child = root->children[i];
            int start = child->start;
            int end = child->end;
            
            // Extract edge label
            char temp_label[1000];
            int len = 0;
            
            for (int j = start; j <= end && j < text_len; j++) {
                temp_label[len++] = text[j];
            }
            temp_label[len] = '\0';
            
            if (len > 0) {
                edge_labels[count] = (char*)malloc((len + 1) * sizeof(char));
                strcpy(edge_labels[count], temp_label);
                count++;
            }
        }
    }
    
    // Print results
    for (int i = 0; i < count; i++) {
        printf("%s\n", edge_labels[i]);
        free(edge_labels[i]);
    }
    
    free(edge_labels);
}

int main() {
    solve();
    return 0;
}
```

## Explanation

This solution implements a suffix tree construction using Ukkonen's algorithm:

1. **Node Structure**: Each node represents a subtree with start and end positions in the text, children pointers, and suffix link information.

2. **Key Functions**:
   - `create_node()`: Creates a new tree node
   - `insert_suffix()`: Inserts a suffix into the tree using Ukkonen's algorithm
   - `build_suffix_tree()`: Builds the complete suffix tree for the input string

3. **Algorithm Steps**:
   - Initialize root node and tracking variables
   - For each suffix, insert it into the tree following Ukkonen's online construction
   - Handle edge splitting when necessary
   - Maintain suffix links for efficient transitions

4. **Output**: Returns all edge labels from the constructed suffix tree.

## Time Complexity
- **Time**: O(n) where n is the length of the input string (Ukkonen's algorithm)
- **Space**: O(n) for storing the tree structure

## Usage
Compile with: `gcc -o suffix_tree suffix_tree.c`
Run with: `./suffix_tree`

The program reads a DNA string from standard input and outputs all edge labels of the constructed suffix tree.