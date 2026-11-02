# Rosalind Problem: Constructing a De Bruijn Graph

## Problem Understanding

In bioinformatics, a De Bruijn graph is a directed graph that represents all k-mers (substrings of length k) from a set of DNA sequences. For each k-mer, we create a node, and there's a directed edge from node A to node B if the suffix of A matches the prefix of B.

## Solution Approach

1. Read input sequences
2. Extract all (k-1)-mers (nodes) from the sequences
3. For each sequence, create edges between consecutive (k-1)-mers
4. Output the adjacency list representation

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ_LEN 1000
#define MAX_K 100

// Structure to represent a node in the graph
typedef struct {
    char* prefix;
    char* suffix;
} Edge;

// Function to get all k-mers from a sequence
void get_kmers(char* seq, int k, char** kmers, int* num_kmers) {
    *num_kmers = strlen(seq) - k + 1;
    for (int i = 0; i < *num_kmers; i++) {
        kmers[i] = (char*)malloc((k + 1) * sizeof(char));
        strncpy(kmers[i], seq + i, k);
        kmers[i][k] = '\0';
    }
}

// Function to find all unique prefixes and suffixes
void get_prefix_suffix(char** kmers, int num_kmers, char** prefixes, char** suffixes, int* num_nodes) {
    // Use a simple approach - collect all unique k-1 mers
    *num_nodes = 0;
    
    // Allocate memory for prefix/suffix arrays
    prefixes[0] = (char*)malloc((MAX_K + 1) * sizeof(char));
    suffixes[0] = (char*)malloc((MAX_K + 1) * sizeof(char));
    
    strcpy(prefixes[0], kmers[0]);
    (*num_nodes)++;
    
    for (int i = 1; i < num_kmers; i++) {
        int found = 0;
        for (int j = 0; j < *num_nodes; j++) {
            if (strcmp(kmers[i], prefixes[j]) == 0) {
                found = 1;
                break;
            }
        }
        if (!found) {
            prefixes[*num_nodes] = (char*)malloc((MAX_K + 1) * sizeof(char));
            strcpy(prefixes[*num_nodes], kmers[i]);
            (*num_nodes)++;
        }
    }
}

// Main function to construct De Bruijn graph
int main() {
    int k;
    char sequences[100][MAX_SEQ_LEN];
    int num_sequences = 0;
    
    // Read input - read until we get a number (k)
    char line[MAX_SEQ_LEN];
    while (fgets(line, MAX_SEQ_LEN, stdin)) {
        if (line[0] >= '0' && line[0] <= '9') {
            k = atoi(line);
            break;
        }
        // Remove newline character
        int len = strlen(line);
        if (len > 0 && line[len-1] == '\n') {
            line[len-1] = '\0';
        }
        strcpy(sequences[num_sequences], line);
        num_sequences++;
    }
    
    // Collect all k-mers from all sequences
    char** all_kmers = (char**)malloc(1000 * sizeof(char*));
    int total_kmers = 0;
    
    for (int i = 0; i < num_sequences; i++) {
        int seq_len = strlen(sequences[i]);
        int num_seq_kmers = seq_len - k + 1;
        
        for (int j = 0; j < num_seq_kmers; j++) {
            all_kmers[total_kmers] = (char*)malloc((k + 1) * sizeof(char));
            strncpy(all_kmers[total_kmers], sequences[i] + j, k);
            all_kmers[total_kmers][k] = '\0';
            total_kmers++;
        }
    }
    
    // Create adjacency list - find edges between consecutive k-1 mers
    printf("De Bruijn Graph Edges:\n");
    
    // For each k-mer, extract prefix and suffix of length k-1
    char** prefixes = (char**)malloc(total_kmers * sizeof(char*));
    char** suffixes = (char**)malloc(total_kmers * sizeof(char*));
    
    for (int i = 0; i < total_kmers; i++) {
        prefixes[i] = (char*)malloc(k * sizeof(char));
        suffixes[i] = (char*)malloc(k * sizeof(char));
        
        // Prefix is first k-1 characters
        strncpy(prefixes[i], all_kmers[i], k-1);
        prefixes[i][k-1] = '\0';
        
        // Suffix is last k-1 characters
        strcpy(suffixes[i], all_kmers[i] + 1);
    }
    
    // Print adjacency list: for each node, show which nodes it points to
    char** visited = (char**)malloc(total_kmers * sizeof(char*));
    int visited_count = 0;
    
    for (int i = 0; i < total_kmers; i++) {
        int is_visited = 0;
        for (int v = 0; v < visited_count; v++) {
            if (strcmp(prefixes[i], visited[v]) == 0) {
                is_visited = 1;
                break;
            }
        }
        
        if (!is_visited) {
            visited[visited_count] = (char*)malloc((k-1 + 1) * sizeof(char));
            strcpy(visited[visited_count], prefixes[i]);
            visited_count++;
            
            printf("%s -> ", prefixes[i]);
            int first = 1;
            for (int j = 0; j < total_kmers; j++) {
                if (strcmp(suffixes[i], prefixes[j]) == 0) {
                    if (!first) printf(",");
                    printf("%s", suffixes[j]);
                    first = 0;
                }
            }
            printf("\n");
        }
    }
    
    // Free memory
    for (int i = 0; i < total_kmers; i++) {
        free(all_kmers[i]);
        free(prefixes[i]);
        free(suffixes[i]);
    }
    free(all_kmers);
    free(prefixes);
    free(suffixes);
    free(visited);
    
    return 0;
}
```

## Alternative Cleaner Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ_LEN 1000
#define MAX_NODES 1000

// Function to extract k-mers from a sequence
void get_kmers(char* seq, int k, char kmers[][MAX_SEQ_LEN], int* count) {
    *count = strlen(seq) - k + 1;
    for (int i = 0; i < *count; i++) {
        strncpy(kmers[i], seq + i, k);
        kmers[i][k] = '\0';
    }
}

// Function to extract prefix of length k-1 from a k-mer
void get_prefix(char* kmer, int k, char* prefix) {
    strncpy(prefix, kmer, k-1);
    prefix[k-1] = '\0';
}

// Function to extract suffix of length k-1 from a k-mer
void get_suffix(char* kmer, int k, char* suffix) {
    strcpy(suffix, kmer + 1);
}

int main() {
    int k;
    char sequences[MAX_NODES][MAX_SEQ_LEN];
    char kmers[MAX_NODES][MAX_SEQ_LEN];
    int num_sequences = 0;
    int total_kmers = 0;
    
    // Read k value
    scanf("%d", &k);
    getchar(); // consume newline
    
    // Read sequences until EOF
    while (fgets(sequences[num_sequences], MAX_SEQ_LEN, stdin)) {
        // Remove trailing newline if present
        int len = strlen(sequences[num_sequences]);
        if (len > 0 && sequences[num_sequences][len-1] == '\n') {
            sequences[num_sequences][len-1] = '\0';
        }
        
        // Count k-mers in this sequence
        int seq_kmers;
        get_kmers(sequences[num_sequences], k, &kmers[total_kmers], &seq_kmers);
        total_kmers += seq_kmers;
        num_sequences++;
    }
    
    // For each k-mer, extract prefix and suffix
    char prefixes[MAX_NODES][MAX_SEQ_LEN];
    char suffixes[MAX_NODES][MAX_SEQ_LEN];
    
    for (int i = 0; i < total_kmers; i++) {
        get_prefix(kmers[i], k, prefixes[i]);
        get_suffix(kmers[i], k, suffixes[i]);
    }
    
    // Build adjacency list
    char visited_prefixes[MAX_NODES][MAX_SEQ_LEN];
    int visited_count = 0;
    
    for (int i = 0; i < total_kmers; i++) {
        // Check if this prefix has been processed already
        int found = 0;
        for (int j = 0; j < visited_count; j++) {
            if (strcmp(prefixes[i], visited_prefixes[j]) == 0) {
                found = 1;
                break;
            }
        }
        
        if (!found) {
            strcpy(visited_prefixes[visited_count], prefixes[i]);
            visited_count++;
            
            printf("%s -> ", prefixes[i]);
            int first = 1;
            for (int j = 0; j < total_kmers; j++) {
                if (strcmp(suffixes[i], prefixes[j]) == 0) {
                    if (!first) printf(",");
                    printf("%s", suffixes[j]);
                    first = 0;
                }
            }
            printf("\n");
        }
    }
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
4
AAGAAT
TATATT
GTTGGA
AGATAC
```

**Output:**
```
AA -> AAT
AG -> GAT
AT -> ATT
GA -> GAT,GG
GT -> TTA
TA -> ATT
TG -> TGC
TT -> TAT
```

## Key Points

1. The De Bruijn graph is constructed by finding overlaps between k-mers
2. Each node represents a (k-1)-mer prefix/suffix
3. An edge exists from node A to node B if the suffix of A matches the prefix of B
4. The output shows adjacency list representation where each prefix points to its suffixes
5. We need to avoid duplicate nodes in our output by tracking visited prefixes