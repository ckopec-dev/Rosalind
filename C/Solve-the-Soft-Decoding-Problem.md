# Rosalind Problem: Solve_the_Soft_Decoding_Problem

## Problem Description
The soft decoding problem involves finding the most likely sequence of hidden states given a sequence of observations and a Hidden Markov Model (HMM). This is typically solved using the Viterbi algorithm.

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_STATES 100
#define MAX_OBSERVATIONS 100
#define MAX_LEN 1000

// Global variables for HMM parameters
double transition_matrix[MAX_STATES][MAX_STATES];
double emission_matrix[MAX_STATES][MAX_STATES];
double initial_probs[MAX_STATES];
int num_states, num_observations;
char observations[MAX_LEN];

// Function to compute the most likely state sequence using Viterbi algorithm
void viterbi() {
    // Viterbi table: viterbi[i][j] = probability of best path ending at state j at time i
    double viterbi[MAX_LEN][MAX_STATES];
    
    // Path tracking for backtracking
    int path[MAX_LEN][MAX_STATES];
    
    int len = strlen(observations);
    
    // Initialize first time step
    for (int i = 0; i < num_states; i++) {
        viterbi[0][i] = initial_probs[i] * emission_matrix[i][observations[0] - 'A'];
        path[0][i] = 0;
    }
    
    // Fill the Viterbi table
    for (int t = 1; t < len; t++) {
        for (int j = 0; j < num_states; j++) {
            double max_prob = -1.0;
            int max_state = 0;
            
            for (int i = 0; i < num_states; i++) {
                double prob = viterbi[t-1][i] * transition_matrix[i][j];
                if (prob > max_prob) {
                    max_prob = prob;
                    max_state = i;
                }
            }
            
            viterbi[t][j] = max_prob * emission_matrix[j][observations[t] - 'A'];
            path[t][j] = max_state;
        }
    }
    
    // Find the most likely final state
    double max_prob = -1.0;
    int best_final_state = 0;
    for (int i = 0; i < num_states; i++) {
        if (viterbi[len-1][i] > max_prob) {
            max_prob = viterbi[len-1][i];
            best_final_state = i;
        }
    }
    
    // Backtrack to find the most likely path
    int result[MAX_LEN];
    result[len-1] = best_final_state;
    for (int t = len - 2; t >= 0; t--) {
        result[t] = path[t+1][result[t+1]];
    }
    
    // Print the result
    for (int i = 0; i < len; i++) {
        printf("%c", result[i] + 'A');
    }
    printf("\n");
}

// Function to read input data
void read_input() {
    // Read number of states and observations
    scanf("%d %d", &num_states, &num_observations);
    
    // Read initial probabilities
    for (int i = 0; i < num_states; i++) {
        scanf("%lf", &initial_probs[i]);
    }
    
    // Read transition matrix
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            scanf("%lf", &transition_matrix[i][j]);
        }
    }
    
    // Read emission matrix
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_observations; j++) {
            scanf("%lf", &emission_matrix[i][j]);
        }
    }
    
    // Read observation sequence
    scanf("%s", observations);
}

int main() {
    read_input();
    viterbi();
    return 0;
}
```

## Explanation

This solution implements the Viterbi algorithm to solve the soft decoding problem:

1. **Input Reading**: The program reads:
   - Number of hidden states and observations
   - Initial state probabilities
   - Transition probability matrix
   - Emission probability matrix
   - Observation sequence

2. **Viterbi Algorithm**:
   - Initialize the first time step with initial probabilities multiplied by emission probabilities
   - For each subsequent time step, calculate the maximum probability of reaching each state from all previous states
   - Track the path that leads to each state for backtracking
   - Backtrack from the most likely final state to reconstruct the optimal path

3. **Output**: The program prints the most likely sequence of hidden states.

## Time Complexity
- Time: O(T × N²) where T is the length of observation sequence and N is the number of states
- Space: O(T × N) for the Viterbi table

## Sample Input Format
```
3 2
0.5 0.3 0.2
0.3 0.5 0.2
0.4 0.1 0.5
0.2 0.6 0.2
0.7 0.3
0.1 0.9
0.8 0.2
ABC
```

The solution handles the standard HMM decoding problem efficiently using dynamic programming principles.