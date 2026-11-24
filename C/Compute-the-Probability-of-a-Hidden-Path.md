# Rosalind Problem: Compute the Probability of a Hidden Path

## Problem Understanding

We need to compute the probability of a hidden path in a Hidden Markov Model (HMM). Given:
- A sequence of observed states (emissions)
- Transition probabilities between hidden states
- Emission probabilities for each state
- Initial state probabilities

The goal is to calculate the probability of the given hidden path.

## Solution Approach

We'll use dynamic programming to compute the probability. For a hidden path of length n, we multiply:
1. Initial state probability
2. Transition probabilities between consecutive states
3. Emission probabilities for each observed state

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_STATES 100
#define MAX_PATH_LENGTH 1000

double compute_hidden_path_probability(char *path, char *states, 
                                       double **transition, double **emission,
                                       double *initial, int num_states) {
    double probability = initial[path[0] - 'A']; // Initial state probability
    
    // Multiply by transition probabilities between consecutive states
    for (int i = 1; i < strlen(path); i++) {
        int prev_state = path[i-1] - 'A';
        int curr_state = path[i] - 'A';
        probability *= transition[prev_state][curr_state];
    }
    
    // Multiply by emission probabilities for each observed state
    for (int i = 0; i < strlen(path); i++) {
        int state = path[i] - 'A';
        int observation = path[i] - 'A'; // Assuming same mapping for observations
        probability *= emission[state][observation];
    }
    
    return probability;
}

int main() {
    // Example usage with sample data
    char path[] = "BABAABAT";
    char states[] = "AB";
    int num_states = 2;
    
    // Initialize matrices
    double **transition = (double**)malloc(num_states * sizeof(double*));
    double **emission = (double**)malloc(num_states * sizeof(double*));
    double *initial = (double*)malloc(num_states * sizeof(double));
    
    for (int i = 0; i < num_states; i++) {
        transition[i] = (double*)malloc(num_states * sizeof(double));
        emission[i] = (double*)malloc(num_states * sizeof(double));
    }
    
    // Sample transition matrix (A->A, A->B, B->A, B->B)
    transition[0][0] = 0.377; transition[0][1] = 0.623;
    transition[1][0] = 0.344; transition[1][1] = 0.656;
    
    // Sample emission matrix (A->A, A->B, B->A, B->B)
    emission[0][0] = 0.592; emission[0][1] = 0.408;
    emission[1][0] = 0.633; emission[1][1] = 0.367;
    
    // Sample initial probabilities
    initial[0] = 0.5; initial[1] = 0.5;
    
    double result = compute_hidden_path_probability(path, states, 
                                                    transition, emission, 
                                                    initial, num_states);
    
    printf("Probability of hidden path: %.10f\n", result);
    
    // Free memory
    for (int i = 0; i < num_states; i++) {
        free(transition[i]);
        free(emission[i]);
    }
    free(transition);
    free(emission);
    free(initial);
    
    return 0;
}
```

## Alternative Implementation with Proper Input Parsing

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double calculate_path_probability(char *path, char **states, 
                                  double **transitions, double **emissions,
                                  double *initial_probs, int num_states) {
    double prob = initial_probs[path[0] - 'A'];
    
    // Multiply by transition probabilities
    for (int i = 1; i < strlen(path); i++) {
        int prev_state = path[i-1] - 'A';
        int curr_state = path[i] - 'A';
        prob *= transitions[prev_state][curr_state];
    }
    
    // Multiply by emission probabilities
    for (int i = 0; i < strlen(path); i++) {
        int state = path[i] - 'A';
        int obs = path[i] - 'A'; // For simplicity, assuming same indexing
        prob *= emissions[state][obs];
    }
    
    return prob;
}

int main() {
    // Read input from stdin or file
    char path[1000];
    char states[100];
    int num_states;
    
    printf("Enter hidden path: ");
    scanf("%s", path);
    
    printf("Enter states (e.g., AB): ");
    scanf("%s", states);
    
    num_states = strlen(states);
    
    // Allocate memory for matrices
    double **transitions = (double**)malloc(num_states * sizeof(double*));
    double **emissions = (double**)malloc(num_states * sizeof(double*));
    double *initial_probs = (double*)malloc(num_states * sizeof(double));
    
    for (int i = 0; i < num_states; i++) {
        transitions[i] = (double*)calloc(num_states, sizeof(double));
        emissions[i] = (double*)calloc(num_states, sizeof(double));
    }
    
    // Read transition matrix
    printf("Enter transition probabilities:\n");
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            scanf("%lf", &transitions[i][j]);
        }
    }
    
    // Read emission matrix
    printf("Enter emission probabilities:\n");
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            scanf("%lf", &emissions[i][j]);
        }
    }
    
    // Read initial probabilities
    printf("Enter initial probabilities:\n");
    for (int i = 0; i < num_states; i++) {
        scanf("%lf", &initial_probs[i]);
    }
    
    double result = calculate_path_probability(path, NULL, 
                                               transitions, emissions, 
                                               initial_probs, num_states);
    
    printf("Probability: %.10f\n", result);
    
    // Free memory
    for (int i = 0; i < num_states; i++) {
        free(transitions[i]);
        free(emissions[i]);
    }
    free(transitions);
    free(emissions);
    free(initial_probs);
    
    return 0;
}
```

## Key Points

1. **Input Format**: The hidden path is given as a string of states
2. **Matrix Dimensions**: Transition and emission matrices are square (n×n)
3. **Indexing**: States are indexed by their position in the alphabet ('A' = 0, 'B' = 1, etc.)
4. **Calculation**: Probability = Initial × ∏(Transition) × ∏(Emission)
5. **Memory Management**: Proper allocation and deallocation of matrices

## Time Complexity
- Time: O(n) where n is the length of the path
- Space: O(n²) for the transition and emission matrices

This solution correctly computes the probability of a given hidden path in a Hidden Markov Model using dynamic programming principles.