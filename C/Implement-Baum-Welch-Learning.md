# Rosalind Problem: Implement Baum-Welch Learning (C Solution)

## Problem Understanding

The Baum-Welch algorithm is used to learn the parameters of a Hidden Markov Model (HMM) from observed sequences when the hidden states are not known. This problem requires implementing the Baum-Welch learning algorithm to estimate transition and emission probabilities.

## Approach

We'll implement the Baum-Welch algorithm using:
1. Forward-backward procedure for computing alpha and beta values
2. Expectation step (computing gamma and xi values)
3. Maximization step (updating transition and emission probabilities)

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_STATES 100
#define MAX_OBSERVATIONS 100
#define MAX_SEQUENCE_LENGTH 1000

// Structure to represent HMM parameters
typedef struct {
    int num_states;
    int num_observations;
    double **transition;    // A[i][j] = P(state j | state i)
    double **emission;      // B[i][k] = P(observation k | state i)
    double *initial;        // pi[i] = P(initial state is i)
} HMM;

// Forward algorithm
void forward(HMM *hmm, char *sequence, int seq_length, double **alpha) {
    // Initialize alpha for first time step
    for (int i = 0; i < hmm->num_states; i++) {
        alpha[0][i] = hmm->initial[i] * hmm->emission[i][sequence[0] - 'A'];
    }
    
    // Forward recursion
    for (int t = 1; t < seq_length; t++) {
        for (int j = 0; j < hmm->num_states; j++) {
            alpha[t][j] = 0.0;
            for (int i = 0; i < hmm->num_states; i++) {
                alpha[t][j] += alpha[t-1][i] * hmm->transition[i][j];
            }
            alpha[t][j] *= hmm->emission[j][sequence[t] - 'A'];
        }
    }
}

// Backward algorithm
void backward(HMM *hmm, char *sequence, int seq_length, double **beta) {
    // Initialize beta for last time step
    for (int i = 0; i < hmm->num_states; i++) {
        beta[seq_length-1][i] = 1.0;
    }
    
    // Backward recursion
    for (int t = seq_length - 2; t >= 0; t--) {
        for (int i = 0; i < hmm->num_states; i++) {
            beta[t][i] = 0.0;
            for (int j = 0; j < hmm->num_states; j++) {
                beta[t][i] += hmm->transition[i][j] * hmm->emission[j][sequence[t+1] - 'A'] * beta[t+1][j];
            }
        }
    }
}

// Baum-Welch learning
void baum_welch(HMM *hmm, char **sequences, int num_sequences, int max_iterations, double epsilon) {
    double **alpha = (double**)malloc(MAX_SEQUENCE_LENGTH * sizeof(double*));
    double **beta = (double**)malloc(MAX_SEQUENCE_LENGTH * sizeof(double*));
    double ***gamma = (double***)malloc(MAX_SEQUENCE_LENGTH * sizeof(double**));
    double ****xi = (double****)malloc(MAX_SEQUENCE_LENGTH * sizeof(double***));
    
    for (int i = 0; i < MAX_SEQUENCE_LENGTH; i++) {
        alpha[i] = (double*)calloc(hmm->num_states, sizeof(double));
        beta[i] = (double*)calloc(hmm->num_states, sizeof(double));
        gamma[i] = (double**)malloc(hmm->num_states * sizeof(double*));
        xi[i] = (double***)malloc(hmm->num_states * sizeof(double**));
        for (int j = 0; j < hmm->num_states; j++) {
            gamma[i][j] = (double*)calloc(hmm->num_states, sizeof(double));
            xi[i][j] = (double**)malloc(hmm->num_states * sizeof(double*));
            for (int k = 0; k < hmm->num_states; k++) {
                xi[i][j][k] = (double*)calloc(1, sizeof(double));
            }
        }
    }
    
    double **new_transition = (double**)malloc(hmm->num_states * sizeof(double*));
    double **new_emission = (double**)malloc(hmm->num_states * sizeof(double*));
    double *new_initial = (double*)calloc(hmm->num_states, sizeof(double));
    
    for (int i = 0; i < hmm->num_states; i++) {
        new_transition[i] = (double*)calloc(hmm->num_states, sizeof(double));
        new_emission[i] = (double*)calloc(hmm->num_observations, sizeof(double));
    }
    
    int iter = 0;
    double old_likelihood = -1e308;
    
    while (iter < max_iterations) {
        double total_likelihood = 0.0;
        
        for (int seq_idx = 0; seq_idx < num_sequences; seq_idx++) {
            char *sequence = sequences[seq_idx];
            int seq_length = strlen(sequence);
            
            // Forward and backward algorithms
            forward(hmm, sequence, seq_length, alpha);
            backward(hmm, sequence, seq_length, beta);
            
            // Compute likelihood for this sequence
            double likelihood = 0.0;
            for (int i = 0; i < hmm->num_states; i++) {
                likelihood += alpha[seq_length-1][i];
            }
            total_likelihood += log(likelihood);
            
            // Compute gamma and xi values
            for (int t = 0; t < seq_length; t++) {
                double norm = 0.0;
                for (int i = 0; i < hmm->num_states; i++) {
                    gamma[t][i] = alpha[t][i] * beta[t][i];
                    norm += gamma[t][i];
                }
                
                // Normalize gamma
                for (int i = 0; i < hmm->num_states; i++) {
                    if (norm > 1e-308) {
                        gamma[t][i] /= norm;
                    } else {
                        gamma[t][i] = 0.0;
                    }
                }
                
                // Compute xi for t < seq_length - 1
                if (t < seq_length - 1) {
                    double norm_xi = 0.0;
                    for (int i = 0; i < hmm->num_states; i++) {
                        for (int j = 0; j < hmm->num_states; j++) {
                            xi[t][i][j] = alpha[t][i] * hmm->transition[i][j] * 
                                         hmm->emission[j][sequence[t+1] - 'A'] * beta[t+1][j];
                            norm_xi += xi[t][i][j];
                        }
                    }
                    
                    // Normalize xi
                    for (int i = 0; i < hmm->num_states; i++) {
                        for (int j = 0; j < hmm->num_states; j++) {
                            if (norm_xi > 1e-308) {
                                xi[t][i][j] /= norm_xi;
                            } else {
                                xi[t][i][j] = 0.0;
                            }
                        }
                    }
                }
            }
            
            // Update initial probabilities
            for (int i = 0; i < hmm->num_states; i++) {
                new_initial[i] += gamma[0][i];
            }
            
            // Update transition probabilities
            for (int i = 0; i < hmm->num_states; i++) {
                double den = 0.0;
                for (int t = 0; t < seq_length - 1; t++) {
                    den += gamma[t][i];
                }
                
                for (int j = 0; j < hmm->num_states; j++) {
                    double num = 0.0;
                    for (int t = 0; t < seq_length - 1; t++) {
                        num += xi[t][i][j];
                    }
                    new_transition[i][j] += num / den;
                }
            }
            
            // Update emission probabilities
            for (int i = 0; i < hmm->num_states; i++) {
                double den = 0.0;
                for (int t = 0; t < seq_length; t++) {
                    den += gamma[t][i];
                }
                
                for (int k = 0; k < hmm->num_observations; k++) {
                    double num = 0.0;
                    for (int t = 0; t < seq_length; t++) {
                        if (sequence[t] - 'A' == k) {
                            num += gamma[t][i];
                        }
                    }
                    new_emission[i][k] += num / den;
                }
            }
        }
        
        // Normalize updated probabilities
        double norm_initial = 0.0;
        for (int i = 0; i < hmm->num_states; i++) {
            norm_initial += new_initial[i];
        }
        
        for (int i = 0; i < hmm->num_states; i++) {
            if (norm_initial > 1e-308) {
                hmm->initial[i] = new_initial[i] / norm_initial;
            } else {
                hmm->initial[i] = 1.0 / hmm->num_states;
            }
        }
        
        for (int i = 0; i < hmm->num_states; i++) {
            double norm_transition = 0.0;
            for (int j = 0; j < hmm->num_states; j++) {
                norm_transition += new_transition[i][j];
            }
            
            for (int j = 0; j < hmm->num_states; j++) {
                if (norm_transition > 1e-308) {
                    hmm->transition[i][j] = new_transition[i][j] / norm_transition;
                } else {
                    hmm->transition[i][j] = 1.0 / hmm->num_states;
                }
            }
        }
        
        for (int i = 0; i < hmm->num_states; i++) {
            double norm_emission = 0.0;
            for (int k = 0; k < hmm->num_observations; k++) {
                norm_emission += new_emission[i][k];
            }
            
            for (int k = 0; k < hmm->num_observations; k++) {
                if (norm_emission > 1e-308) {
                    hmm->emission[i][k] = new_emission[i][k] / norm_emission;
                } else {
                    hmm->emission[i][k] = 1.0 / hmm->num_observations;
                }
            }
        }
        
        // Check convergence
        if (iter > 0) {
            double diff = fabs(total_likelihood - old_likelihood);
            if (diff < epsilon) break;
        }
        old_likelihood = total_likelihood;
        iter++;
    }
    
    // Free memory
    for (int i = 0; i < MAX_SEQUENCE_LENGTH; i++) {
        free(alpha[i]);
        free(beta[i]);
        for (int j = 0; j < hmm->num_states; j++) {
            free(gamma[i][j]);
            for (int k = 0; k < hmm->num_states; k++) {
                free(xi[i][j][k]);
            }
            free(xi[i][j]);
        }
        free(gamma[i]);
        free(xi[i]);
    }
    free(alpha);
    free(beta);
    free(gamma);
    free(xi);
    
    for (int i = 0; i < hmm->num_states; i++) {
        free(new_transition[i]);
        free(new_emission[i]);
    }
    free(new_transition);
    free(new_emission);
    free(new_initial);
}

// Print HMM parameters
void print_hmm(HMM *hmm) {
    printf("Initial probabilities:\n");
    for (int i = 0; i < hmm->num_states; i++) {
        printf("%.3f ", hmm->initial[i]);
    }
    printf("\n\n");
    
    printf("Transition matrix:\n");
    for (int i = 0; i < hmm->num_states; i++) {
        for (int j = 0; j < hmm->num_states; j++) {
            printf("%.3f ", hmm->transition[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("Emission matrix:\n");
    for (int i = 0; i < hmm->num_states; i++) {
        for (int j = 0; j < hmm->num_observations; j++) {
            printf("%.3f ", hmm->emission[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Main function
int main() {
    // Initialize HMM with example parameters
    HMM hmm;
    hmm.num_states = 2;
    hmm.num_observations = 3;
    
    // Allocate memory for HMM parameters
    hmm.transition = (double**)malloc(hmm.num_states * sizeof(double*));
    hmm.emission = (double**)malloc(hmm.num_states * sizeof(double*));
    hmm.initial = (double*)calloc(hmm.num_states, sizeof(double));
    
    for (int i = 0; i < hmm.num_states; i++) {
        hmm.transition[i] = (double*)calloc(hmm.num_states, sizeof(double));
        hmm.emission[i] = (double*)calloc(hmm.num_observations, sizeof(double));
    }
    
    // Initialize with example values
    // Initial probabilities
    hmm.initial[0] = 0.5;
    hmm.initial[1] = 0.5;
    
    // Transition probabilities
    hmm.transition[0][0] = 0.7;
    hmm.transition[0][1] = 0.3;
    hmm.transition[1][0] = 0.4;
    hmm.transition[1][1] = 0.6;
    
    // Emission probabilities
    hmm.emission[0][0] = 0.1;
    hmm.emission[0][1] = 0.4;
    hmm.emission[0][2] = 0.5;
    hmm.emission[1][0] = 0.7;
    hmm.emission[1][1] = 0.2;
    hmm.emission[1][2] = 0.1;
    
    // Example sequences
    char *sequences[] = {"ABC", "BAC"};
    int num_sequences = 2;
    
    printf("Initial HMM parameters:\n");
    print_hmm(&hmm);
    
    // Run Baum-Welch learning
    baum_welch(&hmm, sequences, num_sequences, 100, 1e-6);
    
    printf("Final HMM parameters after Baum-Welch learning:\n");
    print_hmm(&hmm);
    
    // Free memory
    for (int i = 0; i < hmm.num_states; i++) {
        free(hmm.transition[i]);
        free(hmm.emission[i]);
    }
    free(hmm.transition);
    free(hmm.emission);
    free(hmm.initial);
    
    return 0;
}
```

## Key Components of the Solution

1. **Forward Algorithm**: Computes alpha values (probability of being in state i at time t, given observations up to time t)

2. **Backward Algorithm**: Computes beta values (probability of observing future observations from time t+1 onwards, given state i at time t)

3. **Expectation Step**: Uses alpha and beta to compute gamma (posterior probability of being in state i at time t) and xi (posterior probability of transitioning from state i to state j at time t)

4. **Maximization Step**: Updates transition and emission probabilities based on expected counts

5. **Convergence Check**: Stops when likelihood improvement is below threshold

## Time Complexity
- O(N × M² × T) per iteration where N = number of sequences, M = number of states, T = average sequence length
- Overall: O(iterations × N × M² × T)

## Space Complexity
- O(M × T) for storing alpha and beta matrices

This implementation follows the standard Baum-Welch algorithm and should work correctly for the Rosalind problem requirements.