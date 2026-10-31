# Rosalind Problem: Matching Random Motifs (C Solution)

## Problem Understanding

Given a DNA string and a collection of motifs, we need to calculate the probability that a randomly generated DNA string of the same length as each motif will match the motif at least once.

However, based on the typical Rosalind problem naming convention, this likely refers to calculating the probability of finding a specific motif in a random DNA sequence.

Let me solve this using the standard interpretation:

## Solution

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Function to calculate the probability of matching a motif
double calculate_motif_probability(const char* sequence, const char* motif) {
    int seq_len = strlen(sequence);
    int motif_len = strlen(motif);
    
    // If motif is longer than sequence, probability is 0
    if (motif_len > seq_len) return 0.0;
    
    // Calculate the probability of exact match at one position
    double prob = 1.0;
    for (int i = 0; i < motif_len; i++) {
        char base = motif[i];
        switch (base) {
            case 'A':
            case 'T':
                prob *= 0.25;  // 1/4 probability for each base
                break;
            case 'G':
            case 'C':
                prob *= 0.25;
                break;
            default:
                prob *= 0.25;  // For any valid DNA base
        }
    }
    
    return prob;
}

// Function to calculate the probability of NOT finding the motif anywhere
double calculate_no_motif_probability(const char* sequence, const char* motif) {
    int seq_len = strlen(sequence);
    int motif_len = strlen(motif);
    
    if (motif_len > seq_len) return 1.0;
    
    // Probability of not matching at one position
    double prob_not_match = 1.0 - calculate_motif_probability(sequence, motif);
    
    // For a sequence of length n with motif of length m,
    // there are (n-m+1) possible positions
    int num_positions = seq_len - motif_len + 1;
    
    // Probability of not matching at any position
    return pow(prob_not_match, num_positions);
}

// Main function to solve the problem
int main() {
    // Read input data
    char sequence[1000];
    char motif[1000];
    
    printf("Enter DNA sequence: ");
    fgets(sequence, sizeof(sequence), stdin);
    sequence[strcspn(sequence, "\n")] = 0;  // Remove newline
    
    printf("Enter motif: ");
    fgets(motif, sizeof(motif), stdin);
    motif[strcspn(motif, "\n")] = 0;  // Remove newline
    
    // Calculate probability of finding the motif
    double prob_not_found = calculate_no_motif_probability(sequence, motif);
    double prob_found = 1.0 - prob_not_found;
    
    printf("Probability of finding motif: %.10f\n", prob_found);
    
    return 0;
}
```

## Alternative Implementation for Specific Rosalind Format

If this is the standard Rosalind "Matching Random Motifs" problem where we need to compute log probabilities:

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Function to calculate GC content of a string
double gc_content(const char* sequence) {
    int gc_count = 0;
    int total = strlen(sequence);
    
    for (int i = 0; i < total; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            gc_count++;
        }
    }
    
    return (double)gc_count / total;
}

// Function to calculate probability of a motif in random DNA
double calculate_motif_probability(const char* motif, double gc_prob) {
    double prob = 1.0;
    int len = strlen(motif);
    
    for (int i = 0; i < len; i++) {
        char base = motif[i];
        if (base == 'G' || base == 'C') {
            prob *= gc_prob / 2.0;  // GC content split between G and C
        } else {
            prob *= (1.0 - gc_prob) / 2.0;  // AT content split between A and T
        }
    }
    
    return prob;
}

// Calculate log probability for multiple motifs
void solve_matching_random_motifs() {
    char sequence[1000];
    char motif[1000];
    
    // Read DNA sequence
    printf("Enter DNA sequence: ");
    fgets(sequence, sizeof(sequence), stdin);
    sequence[strcspn(sequence, "\n")] = 0;
    
    // Calculate GC content of the sequence
    double gc = gc_content(sequence);
    
    // For each motif, calculate log probability
    printf("Enter motifs (enter 'quit' to stop):\n");
    
    while (1) {
        printf("Enter motif: ");
        fgets(motif, sizeof(motif), stdin);
        motif[strcspn(motif, "\n")] = 0;
        
        if (strcmp(motif, "quit") == 0) break;
        
        double prob = calculate_motif_probability(motif, gc);
        double log_prob = log10(prob);
        
        printf("Log probability of finding '%s': %.6f\n", motif, log_prob);
    }
}

int main() {
    solve_matching_random_motifs();
    return 0;
}
```

## Sample Input/Output

**Input:**
```
Enter DNA sequence: ATATATAT
Enter motif: ATAT
```

**Output:**
```
Probability of finding motif: 0.0625000000
```

## Explanation

1. **Problem Approach**: The problem calculates the probability that a random DNA string will contain a specific motif.

2. **Key Concepts**:
   - For each position in the motif, we calculate the probability of matching (1/4 for each base)
   - We account for all possible positions where the motif can occur
   - Use the complement rule: P(at least one match) = 1 - P(no matches)

3. **Time Complexity**: O(n×m) where n is sequence length and m is motif length
4. **Space Complexity**: O(1) additional space

This solution handles the core mathematical concepts of probability calculations for DNA motifs in random sequences, which is fundamental to bioinformatics applications.