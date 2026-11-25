# Rosalind Problem: Perform a Multiple Sequence Alignment with a Profile HMM

## Problem Understanding

In this problem, we need to implement a multiple sequence alignment using a Profile Hidden Markov Model (Profile HMM). The key components are:
- Input: A set of sequences aligned to a consensus sequence
- Output: A profile HMM that represents the alignment
- We'll use the Viterbi algorithm to find the optimal path through the HMM

## Solution Approach

I'll implement a Profile HMM with:
1. Three main states (M, D, I) for each position in the consensus
2. Transition probabilities between states
3. Emission probabilities for each symbol
4. Viterbi algorithm to find optimal alignment

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_SEQ_LEN 1000
#define MAX_NUM_SEQS 100
#define MAX_STATES 300
#define ALPHABET_SIZE 26

// State types
typedef enum {
    MATCH = 0,
    DELETE = 1,
    INSERT = 2
} StateType;

// HMM state structure
typedef struct {
    int state_id;
    StateType type;
    int position; // Position in consensus sequence
} HMMState;

// Profile HMM structure
typedef struct {
    int num_states;
    int num_positions;
    double **transition_probs;
    double **emission_probs;
    HMMState *states;
    char *consensus;
    int *seq_lengths;
    int num_sequences;
} ProfileHMM;

// Global variables for sequences
char sequences[MAX_NUM_SEQS][MAX_SEQ_LEN];
int num_sequences = 0;
char consensus_seq[MAX_SEQ_LEN];

// Function prototypes
ProfileHMM* create_profile_hmm(char **seqs, int n_seqs, char *consensus);
void initialize_transition_probs(ProfileHMM *hmm);
void initialize_emission_probs(ProfileHMM *hmm, char **seqs, int n_seqs);
double viterbi_algorithm(ProfileHMM *hmm, char *query_seq, double ***viterbi_matrix);
void free_profile_hmm(ProfileHMM *hmm);

// Create a profile HMM from aligned sequences
ProfileHMM* create_profile_hmm(char **seqs, int n_seqs, char *consensus) {
    ProfileHMM *hmm = (ProfileHMM*)malloc(sizeof(ProfileHMM));
    
    // Initialize basic parameters
    hmm->num_sequences = n_seqs;
    hmm->num_positions = strlen(consensus);
    hmm->num_states = 3 * hmm->num_positions + 3; // M, D, I states plus start/end
    
    // Allocate memory for states
    hmm->states = (HMMState*)malloc(hmm->num_states * sizeof(HMMState));
    
    // Create states (M0, D0, I0, M1, D1, I1, ..., Mn, Dn, In)
    int state_id = 0;
    for (int i = 0; i < hmm->num_positions; i++) {
        hmm->states[state_id].state_id = state_id;
        hmm->states[state_id].type = MATCH;
        hmm->states[state_id].position = i;
        state_id++;
        
        hmm->states[state_id].state_id = state_id;
        hmm->states[state_id].type = DELETE;
        hmm->states[state_id].position = i;
        state_id++;
        
        hmm->states[state_id].state_id = state_id;
        hmm->states[state_id].type = INSERT;
        hmm->states[state_id].position = i;
        state_id++;
    }
    
    // Add final states
    hmm->states[state_id].state_id = state_id;
    hmm->states[state_id].type = MATCH;
    hmm->states[state_id].position = -1; // End state
    
    // Allocate transition probabilities
    hmm->transition_probs = (double**)malloc(hmm->num_states * sizeof(double*));
    for (int i = 0; i < hmm->num_states; i++) {
        hmm->transition_probs[i] = (double*)calloc(hmm->num_states, sizeof(double));
    }
    
    // Allocate emission probabilities
    hmm->emission_probs = (double**)malloc(hmm->num_states * sizeof(double*));
    for (int i = 0; i < hmm->num_states; i++) {
        hmm->emission_probs[i] = (double*)calloc(ALPHABET_SIZE, sizeof(double));
    }
    
    // Initialize probabilities
    initialize_transition_probs(hmm);
    initialize_emission_probs(hmm, seqs, n_seqs);
    
    return hmm;
}

// Initialize transition probabilities
void initialize_transition_probs(ProfileHMM *hmm) {
    // Set up basic transitions between consecutive states
    for (int i = 0; i < hmm->num_positions; i++) {
        int m_idx = 3 * i;
        int d_idx = 3 * i + 1;
        int i_idx = 3 * i + 2;
        
        // From M_i to M_{i+1} (match)
        if (i < hmm->num_positions - 1) {
            hmm->transition_probs[m_idx][m_idx + 3] = 0.7;  // Match to match
            hmm->transition_probs[m_idx][d_idx + 3] = 0.1;  // Match to delete
            hmm->transition_probs[m_idx][i_idx + 3] = 0.2;  // Match to insert
        } else {
            hmm->transition_probs[m_idx][hmm->num_states - 1] = 1.0; // End state
        }
        
        // From D_i to M_{i+1} (delete to match)
        if (i < hmm->num_positions - 1) {
            hmm->transition_probs[d_idx][m_idx + 3] = 0.9;
            hmm->transition_probs[d_idx][d_idx + 3] = 0.1;
        } else {
            hmm->transition_probs[d_idx][hmm->num_states - 1] = 1.0; // End state
        }
        
        // From I_i to M_{i+1} (insert to match)
        if (i < hmm->num_positions - 1) {
            hmm->transition_probs[i_idx][m_idx + 3] = 0.9;
            hmm->transition_probs[i_idx][i_idx + 3] = 0.1;
        } else {
            hmm->transition_probs[i_idx][hmm->num_states - 1] = 1.0; // End state
        }
    }
}

// Initialize emission probabilities based on observed sequences
void initialize_emission_probs(ProfileHMM *hmm, char **seqs, int n_seqs) {
    // For simplicity, we'll use a basic approach:
    // Count how often each symbol appears in each position of the consensus
    
    // Initialize all emissions to small values for smoothing
    for (int i = 0; i < hmm->num_states; i++) {
        for (int j = 0; j < ALPHABET_SIZE; j++) {
            hmm->emission_probs[i][j] = 1e-6;
        }
    }
    
    // Count emissions from sequences
    for (int s = 0; s < n_seqs; s++) {
        int seq_len = strlen(seqs[s]);
        int pos = 0;
        
        for (int i = 0; i < seq_len && pos < hmm->num_positions; i++) {
            char c = seqs[s][i];
            if (c != '-') { // Skip gaps
                int state_idx = 3 * pos + 0; // Match state
                if (c >= 'A' && c <= 'Z') {
                    hmm->emission_probs[state_idx][c - 'A'] += 1.0;
                }
                pos++;
            }
        }
    }
    
    // Normalize emission probabilities
    for (int i = 0; i < hmm->num_states; i++) {
        double sum = 0.0;
        for (int j = 0; j < ALPHABET_SIZE; j++) {
            sum += hmm->emission_probs[i][j];
        }
        if (sum > 0) {
            for (int j = 0; j < ALPHABET_SIZE; j++) {
                hmm->emission_probs[i][j] /= sum;
            }
        }
    }
}

// Viterbi algorithm implementation
double viterbi_algorithm(ProfileHMM *hmm, char *query_seq, double ***viterbi_matrix) {
    int seq_len = strlen(query_seq);
    int num_states = hmm->num_states;
    
    // Allocate viterbi matrix
    *viterbi_matrix = (double**)malloc(seq_len * sizeof(double*));
    for (int i = 0; i < seq_len; i++) {
        (*viterbi_matrix)[i] = (double*)calloc(num_states, sizeof(double));
    }
    
    // Initialize first row (base case)
    int start_state = 0;
    for (int i = 0; i < num_states; i++) {
        if (i == 0) {
            (*viterbi_matrix)[0][i] = log(1.0 / num_states); // Equal probability for start
        } else {
            (*viterbi_matrix)[0][i] = -INFINITY;
        }
    }
    
    // Fill the viterbi matrix using dynamic programming
    for (int t = 0; t < seq_len; t++) {
        char c = query_seq[t];
        int c_idx = (c >= 'A' && c <= 'Z') ? (c - 'A') : 0;
        
        for (int j = 0; j < num_states; j++) {
            if (t == 0) {
                // For first position, set initial probabilities
                (*viterbi_matrix)[t][j] = log(1.0 / num_states) + 
                                          log(hmm->emission_probs[j][c_idx]);
            } else {
                double max_prob = -INFINITY;
                
                // Find maximum probability from previous states
                for (int i = 0; i < num_states; i++) {
                    if ((*viterbi_matrix)[t-1][i] > -INFINITY) {
                        double prob = (*viterbi_matrix)[t-1][i] + 
                                     log(hmm->transition_probs[i][j]) +
                                     log(hmm->emission_probs[j][c_idx]);
                        if (prob > max_prob) {
                            max_prob = prob;
                        }
                    }
                }
                
                (*viterbi_matrix)[t][j] = max_prob;
            }
        }
    }
    
    // Return the maximum probability to end state
    double result = -INFINITY;
    for (int i = 0; i < num_states; i++) {
        if ((*viterbi_matrix)[seq_len-1][i] > result) {
            result = (*viterbi_matrix)[seq_len-1][i];
        }
    }
    
    return result;
}

// Free memory allocated for Profile HMM
void free_profile_hmm(ProfileHMM *hmm) {
    if (hmm == NULL) return;
    
    // Free transition probabilities
    for (int i = 0; i < hmm->num_states; i++) {
        free(hmm->transition_probs[i]);
    }
    free(hmm->transition_probs);
    
    // Free emission probabilities
    for (int i = 0; i < hmm->num_states; i++) {
        free(hmm->emission_probs[i]);
    }
    free(hmm->emission_probs);
    
    // Free states
    free(hmm->states);
    
    free(hmm);
}

// Main function to solve the problem
int main() {
    // Read input (example format)
    char *seqs[] = {"ACGT", "ACG-", "AC-T", "A--T"};
    int n_seqs = 4;
    char consensus[] = "ACGT";
    
    // Create profile HMM
    ProfileHMM *hmm = create_profile_hmm(seqs, n_seqs, consensus);
    
    // Test with a query sequence
    char query_seq[] = "ACGT";
    double ***viterbi_matrix = NULL;
    
    double score = viterbi_algorithm(hmm, query_seq, &viterbi_matrix);
    
    printf("Viterbi score: %f\n", score);
    
    // Print some information about the HMM
    printf("Profile HMM created with %d states and %d positions\n", 
           hmm->num_states, hmm->num_positions);
    
    // Clean up
    for (int i = 0; i < strlen(query_seq); i++) {
        free(viterbi_matrix[i]);
    }
    free(viterbi_matrix);
    free_profile_hmm(hmm);
    
    return 0;
}
```

## Key Features of the Solution

1. **Profile HMM Structure**: Implements a standard Profile HMM with Match (M), Delete (D), and Insert (I) states for each position.

2. **Transition Probabilities**: Sets up basic transitions between consecutive states with appropriate probabilities.

3. **Emission Probabilities**: Calculates emission probabilities from the aligned sequences.

4. **Viterbi Algorithm**: Implements the dynamic programming algorithm to find the optimal path through the HMM.

5. **Memory Management**: Properly allocates and frees memory for all data structures.

## How to Use

1. Compile with: `gcc -o profile_hmm profile_hmm.c -lm`
2. Modify the input sequences in the main function
3. The program will output the Viterbi score for the alignment

This implementation provides a solid foundation for performing multiple sequence alignment using Profile HMMs, which is essential for bioinformatics applications like protein family classification and structure prediction.