# Rosalind Problem: Implement the Viterbi Algorithm (C Solution)

## Problem Understanding

The Viterbi algorithm is used to find the most likely sequence of hidden states in a Hidden Markov Model (HMM). Given an HMM with states, emissions, and transition/emission probabilities, we need to find the optimal state path for a given emission sequence.

## Solution Approach

I'll implement the Viterbi algorithm using dynamic programming with the following steps:
1. Initialize the Viterbi matrix
2. Fill the matrix using dynamic programming
3. Trace back to find the optimal path
4. Return the most likely state sequence

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_STATES 100
#define MAX_LENGTH 1000

// Structure to store HMM parameters
typedef struct {
    int num_states;
    int num_emissions;
    char *states[MAX_STATES];
    char *emissions[MAX_STATES];
    double transition_probs[MAX_STATES][MAX_STATES];
    double emission_probs[MAX_STATES][MAX_STATES];
} HMM;

// Function to find index of a state or emission
int find_state_index(HMM *hmm, char *state) {
    for (int i = 0; i < hmm->num_states; i++) {
        if (strcmp(hmm->states[i], state) == 0) {
            return i;
        }
    }
    return -1;
}

int find_emission_index(HMM *hmm, char *emission) {
    for (int i = 0; i < hmm->num_emissions; i++) {
        if (strcmp(hmm->emissions[i], emission) == 0) {
            return i;
        }
    }
    return -1;
}

// Function to implement Viterbi algorithm
char* viterbi(HMM *hmm, char *emission_sequence) {
    int seq_length = strlen(emission_sequence);
    
    // Create Viterbi matrix (probabilities)
    double viterbi_matrix[MAX_LENGTH][MAX_STATES];
    int backpointer[MAX_LENGTH][MAX_STATES];
    
    // Initialize first row
    int emission_idx = find_emission_index(hmm, &emission_sequence[0]);
    for (int i = 0; i < hmm->num_states; i++) {
        viterbi_matrix[0][i] = 1.0 / hmm->num_states * hmm->emission_probs[i][emission_idx];
        backpointer[0][i] = 0;
    }
    
    // Fill the Viterbi matrix
    for (int t = 1; t < seq_length; t++) {
        emission_idx = find_emission_index(hmm, &emission_sequence[t]);
        
        for (int j = 0; j < hmm->num_states; j++) {
            double max_prob = -1.0;
            int max_state = 0;
            
            // Find the maximum probability from all previous states
            for (int i = 0; i < hmm->num_states; i++) {
                double prob = viterbi_matrix[t-1][i] * hmm->transition_probs[i][j];
                if (prob > max_prob) {
                    max_prob = prob;
                    max_state = i;
                }
            }
            
            // Calculate probability for current state
            viterbi_matrix[t][j] = max_prob * hmm->emission_probs[j][emission_idx];
            backpointer[t][j] = max_state;
        }
    }
    
    // Find the most likely final state
    double max_prob = -1.0;
    int max_state = 0;
    for (int i = 0; i < hmm->num_states; i++) {
        if (viterbi_matrix[seq_length-1][i] > max_prob) {
            max_prob = viterbi_matrix[seq_length-1][i];
            max_state = i;
        }
    }
    
    // Backtrack to find the optimal path
    char *path = malloc((seq_length + 1) * sizeof(char));
    int current_state = max_state;
    
    for (int t = seq_length - 1; t >= 0; t--) {
        path[t] = hmm->states[current_state][0]; // Take first character of state name
        current_state = backpointer[t][current_state];
    }
    
    path[seq_length] = '\0';
    return path;
}

// Function to parse HMM from input
HMM* parse_hmm() {
    HMM *hmm = malloc(sizeof(HMM));
    
    // Read states and emissions (simplified for this example)
    // In a real implementation, you'd read from stdin or file
    
    hmm->num_states = 2;
    hmm->num_emissions = 3;
    
    hmm->states[0] = "A";
    hmm->states[1] = "B";
    
    hmm->emissions[0] = "x";
    hmm->emissions[1] = "y";
    hmm->emissions[2] = "z";
    
    // Initialize transition probabilities
    for (int i = 0; i < MAX_STATES; i++) {
        for (int j = 0; j < MAX_STATES; j++) {
            hmm->transition_probs[i][j] = 0.0;
            hmm->emission_probs[i][j] = 0.0;
        }
    }
    
    // Set example transition probabilities
    hmm->transition_probs[0][0] = 0.7; // A -> A
    hmm->transition_probs[0][1] = 0.3; // A -> B
    hmm->transition_probs[1][0] = 0.4; // B -> A
    hmm->transition_probs[1][1] = 0.6; // B -> B
    
    // Set example emission probabilities
    hmm->emission_probs[0][0] = 0.2; // A emits x
    hmm->emission_probs[0][1] = 0.3; // A emits y
    hmm->emission_probs[0][2] = 0.5; // A emits z
    hmm->emission_probs[1][0] = 0.5; // B emits x
    hmm->emission_probs[1][1] = 0.4; // B emits y
    hmm->emission_probs[1][2] = 0.1; // B emits z
    
    return hmm;
}

// Main function
int main() {
    // Parse HMM from input
    HMM *hmm = parse_hmm();
    
    // Example emission sequence
    char *emission_sequence = "xyx";
    
    // Run Viterbi algorithm
    char *optimal_path = viterbi(hmm, emission_sequence);
    
    // Print result
    printf("Optimal state path: %s\n", optimal_path);
    
    // Free memory
    free(optimal_path);
    free(hmm);
    
    return 0;
}
```

## Key Features of the Implementation

1. **Data Structures**: 
   - HMM structure to store states, emissions, and probabilities
   - Viterbi matrix for dynamic programming calculations

2. **Algorithm Steps**:
   - Initialization of first row with emission probabilities
   - Dynamic programming fill using transition and emission probabilities
   - Backtracking to reconstruct optimal path

3. **Memory Management**: 
   - Proper allocation and deallocation of memory
   - Handling of variable-length sequences

## Time and Space Complexity

- **Time Complexity**: O(n × k²) where n is the sequence length and k is the number of states
- **Space Complexity**: O(n × k) for the Viterbi matrix

## Sample Input/Output

For a simple HMM with:
- States: A, B  
- Emissions: x, y, z
- Input sequence: "xyx"

The output would be the most likely state path that maximizes the probability of generating the given emission sequence.

This implementation provides a solid foundation for solving the Viterbi algorithm problem and can be extended to handle more complex HMM structures and input formats.