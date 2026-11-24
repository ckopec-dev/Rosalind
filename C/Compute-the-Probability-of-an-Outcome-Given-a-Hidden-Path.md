# Rosalind Problem: Compute the Probability of an Outcome Given a Hidden Path

## Problem Description
Given a hidden path and an emission matrix, compute the probability of observing a specific sequence of outcomes.

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    // Read input data
    char path[1000];
    char states[1000];
    char observations[1000];
    
    // Read the hidden path
    scanf("%s", path);
    
    // Read the states
    scanf("%s", states);
    
    // Read the observations
    scanf("%s", observations);
    
    int path_len = strlen(path);
    int states_len = strlen(states);
    int obs_len = strlen(observations);
    
    // Create emission matrix (assuming 2x2 for simplicity)
    // In practice, this would be read from input
    double emission[10][10];
    
    // Initialize emission probabilities (example values)
    // This should be replaced with actual input reading
    emission[0][0] = 0.5; emission[0][1] = 0.5;
    emission[1][0] = 0.3; emission[1][1] = 0.7;
    
    double probability = 1.0;
    
    // Calculate probability of observation given path
    for (int i = 0; i < obs_len; i++) {
        char current_state = path[i];
        char current_observation = observations[i];
        
        // Convert characters to indices (assuming A=0, B=1, etc.)
        int state_idx = current_state - 'A';
        int obs_idx = current_observation - 'A';
        
        // Multiply by emission probability
        probability *= emission[state_idx][obs_idx];
    }
    
    printf("%.10f\n", probability);
    
    return 0;
}
```

## Alternative Implementation with Proper Input Reading

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    // Read the hidden path
    char path[1000];
    scanf("%s", path);
    
    // Read the states (hidden states)
    char states[1000];
    scanf("%s", states);
    
    // Read the observations
    char observations[1000];
    scanf("%s", observations);
    
    int path_len = strlen(path);
    int obs_len = strlen(observations);
    
    // Read emission matrix dimensions
    int num_states, num_observations;
    scanf("%d %d", &num_states, &num_observations);
    
    // Create emission probability matrix
    double **emission = (double**)malloc(num_states * sizeof(double*));
    for (int i = 0; i < num_states; i++) {
        emission[i] = (double*)malloc(num_observations * sizeof(double));
    }
    
    // Read emission probabilities
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_observations; j++) {
            scanf("%lf", &emission[i][j]);
        }
    }
    
    double probability = 1.0;
    
    // Calculate the probability of the observation sequence given the path
    for (int i = 0; i < obs_len; i++) {
        char current_state = path[i];
        char current_observation = observations[i];
        
        // Convert to indices
        int state_idx = current_state - 'A';
        int obs_idx = current_observation - 'A';
        
        // Multiply by emission probability
        probability *= emission[state_idx][obs_idx];
    }
    
    printf("%.10f\n", probability);
    
    // Free memory
    for (int i = 0; i < num_states; i++) {
        free(emission[i]);
    }
    free(emission);
    
    return 0;
}
```

## Key Points

1. **Input Format**: The program reads a hidden path, states, and observations from standard input
2. **Emission Matrix**: Probability of observing each character given a hidden state
3. **Calculation**: Multiply emission probabilities for each position in the path and observation sequence
4. **Indexing**: Characters are converted to indices (A=0, B=1, etc.) for matrix access

## Time Complexity
O(n) where n is the length of the observation sequence

## Space Complexity  
O(m × k) where m is the number of states and k is the number of observations

This solution computes the probability of observing a specific sequence given a known hidden path using the emission probabilities from the HMM model.