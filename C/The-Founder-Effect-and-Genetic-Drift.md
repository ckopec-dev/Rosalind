# Rosalind Problem: The_Founder_Effect_and_Genetic_Drift

## Problem Understanding

This problem asks us to simulate genetic drift in a population that starts with a founder effect - specifically, we start with 2N individuals who are all heterozygous (Aa) at a single locus. We need to calculate the probability that the allele A will eventually be lost from the population.

## Approach

This is a classic problem in population genetics involving the Wright-Fisher model of genetic drift. For a population of size 2N with initial frequency p of allele A, we want to find the probability that allele A gets lost (frequency becomes 0).

The key insight is that this follows a diffusion process where:
- Starting with 2N individuals, each heterozygous (Aa)
- So initial frequency p = 0.5
- We need to calculate the probability of ultimate fixation or loss

For a Wright-Fisher model starting with k copies of an allele in a population of size 2N, the probability of eventual fixation is k/(2N).

Since we start with N copies of A (from 2N heterozygotes) and N copies of a, the probability that A gets lost is:
P(A lost) = 1 - P(A fixed) = 1 - N/(2N) = 1/2

Wait, let me reconsider. Let's think step by step.

## Solution Approach

Let's use the mathematical solution for genetic drift in a Wright-Fisher model:

If we start with k copies of allele A in a population of 2N individuals, and we want to find the probability that A gets lost (i.e., goes to frequency 0), we can use the formula derived from diffusion theory.

The probability of ultimate loss is:
P(lose) = 1 - (k/(2N))

But wait - let me re-read the problem. It's asking for the probability that the allele A will be lost from the population.

If we start with 2N individuals who are all heterozygous (Aa), then:
- We have N copies of A and N copies of a
- So k = N, 2N = 2N

The probability of A being lost is: 1 - (N/(2N)) = 1/2

Actually, let me think about this more carefully.

## Mathematical Solution

For the Wright-Fisher model with 2N diploid individuals:
- Initial state: all heterozygous (Aa)
- So we start with N copies of A and N copies of a
- Probability that A gets lost = 1 - probability that A gets fixed

In the case where the initial frequency is 0.5, the probability of fixation for either allele is:
P(A fixed) = 0.5 (by symmetry since initial frequency is 0.5)

Therefore, P(A lost) = 1 - 0.5 = 0.5

But let me implement a simulation approach to verify and match the expected output format.

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Function to simulate genetic drift
double simulate_genetic_drift(int N, int num_simulations) {
    int total_individuals = 2 * N;
    double total_loss = 0.0;
    
    for (int sim = 0; sim < num_simulations; sim++) {
        // Start with N copies of A and N copies of a
        int count_A = N;
        int count_a = N;
        
        // Simulate until one allele is lost
        while (count_A > 0 && count_a > 0) {
            // Create new generation by sampling
            int new_count_A = 0;
            
            for (int i = 0; i < total_individuals; i++) {
                // Sample two alleles from current population
                int allele1 = (rand() % (count_A + count_a) < count_A) ? 1 : 0;
                int allele2 = (rand() % (count_A + count_a) < count_A) ? 1 : 0;
                
                // Add to new generation
                if (allele1 == 1 || allele2 == 1) {
                    new_count_A++;
                }
            }
            
            // Update counts for next generation
            count_A = new_count_A;
            count_a = total_individuals - new_count_A;
        }
        
        // If we ended with count_a = 0, A was lost
        if (count_a == 0) {
            total_loss += 1.0;
        }
    }
    
    return total_loss / num_simulations;
}

int main() {
    int N;
    printf("Enter population size parameter N: ");
    scanf("%d", &N);
    
    // Set seed for reproducibility
    srand(42);
    
    // For analytical solution, we know that with initial frequency 0.5,
    // probability of A being lost = 0.5
    double result = 0.5;
    
    printf("Probability that allele A will be lost: %.10f\n", result);
    
    return 0;
}
```

## Alternative Mathematical Approach

Looking at this more carefully, let's consider the exact mathematical solution:

For a population of size 2N with initial frequency p of allele A, and we want probability of loss of A.

When p = 0.5 (which is our case), the probability that A gets lost is:
P(A lost) = 1 - p = 1 - 0.5 = 0.5

This can be proven using the theory of random walks or diffusion equations in population genetics.

```c
#include <stdio.h>

int main() {
    int N;
    printf("Enter population size parameter N: ");
    scanf("%d", &N);
    
    // For a population starting with all heterozygotes (p=0.5),
    // the probability of losing allele A is 0.5
    double probability = 0.5;
    
    printf("%.10f\n", probability);
    
    return 0;
}
```

## Final Answer

The probability that allele A will be lost from the population, starting with 2N heterozygous individuals, is **0.5**.

```c
#include <stdio.h>

int main() {
    int N;
    scanf("%d", &N);
    printf("0.5000000000\n");
    return 0;
}
```

This result makes intuitive sense: since we start with equal frequencies of both alleles (0.5 each), by symmetry, the probability that either allele gets lost should be 0.5.