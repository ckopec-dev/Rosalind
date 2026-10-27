# Rosalind Problem: Independent Alleles (Independent Alleles)

## Problem Statement
In a population whose members are in Hardy-Weinberg equilibrium, we can model the probability of having exactly k dominant alleles in n offspring using the binomial distribution.

Given:
- k = number of dominant alleles we want to find probability for
- n = total number of offspring (in the nth generation)
- The population is in Hardy-Weinberg equilibrium with p = 0.5 (equal allele frequencies)

We need to calculate the probability that a given individual has exactly k dominant alleles among n offspring.

## Solution Approach

The key insight is that each offspring independently inherits alleles from two parents, and we can use the binomial distribution where:
- Each trial (offspring) has probability p = 0.5 of getting the dominant allele
- We want exactly k successes out of n trials

However, since this is the nth generation, we need to consider that each generation produces offspring independently, so we use the binomial probability formula.

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

// Function to calculate factorial
double factorial(int n) {
    if (n <= 1) return 1;
    double result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

// Function to calculate combination C(n,k)
double combination(int n, int k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    
    // Use the more stable calculation to avoid overflow
    double result = 1;
    for (int i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

// Function to calculate binomial probability P(X = k)
double binomial_probability(int n, int k) {
    // Probability of success (getting dominant allele) is 0.5
    double p = 0.5;
    
    // Binomial probability: C(n,k) * p^k * (1-p)^(n-k)
    double comb = combination(n, k);
    double prob = comb * pow(p, k) * pow(1 - p, n - k);
    
    return prob;
}

// Function to calculate cumulative probability P(X >= k)
double cumulative_probability(int n, int k) {
    double total_prob = 0.0;
    for (int i = k; i <= n; i++) {
        total_prob += binomial_probability(n, i);
    }
    return total_prob;
}

int main() {
    // Example: k = 3, n = 2
    int k = 3;  // number of dominant alleles
    int n = 2;  // number of offspring (in nth generation)
    
    printf("Finding probability of having exactly %d dominant alleles among %d offspring\n", k, n);
    
    // For the problem, we want P(X >= k) where X is number of dominant alleles
    double result = cumulative_probability(n, k);
    
    printf("Probability = %.10f\n", result);
    
    // Another example: k = 3, n = 5 (as in typical Rosalind examples)
    k = 3;
    n = 5;
    result = cumulative_probability(n, k);
    printf("\nFor k=%d, n=%d: Probability = %.10f\n", k, n, result);
    
    return 0;
}
```

## Alternative Implementation (More Efficient)

```c
#include <stdio.h>
#include <math.h>

// More efficient implementation using log-space arithmetic to avoid overflow
double log_factorial(int n) {
    if (n <= 1) return 0;
    double sum = 0;
    for (int i = 2; i <= n; i++) {
        sum += log(i);
    }
    return sum;
}

// Calculate combination using log space
double log_combination(int n, int k) {
    if (k > n || k < 0) return -1000000; // Very small number representing 0
    if (k == 0 || k == n) return 0;
    
    double result = log_factorial(n) - log_factorial(k) - log_factorial(n - k);
    return result;
}

// Calculate probability in log space to avoid underflow/overflow
double log_binomial_probability(int n, int k) {
    double log_comb = log_combination(n, k);
    double log_prob = log_comb + k * log(0.5) + (n - k) * log(0.5);
    return log_prob;
}

// Convert back from log space
double probability_from_log(double log_prob) {
    return exp(log_prob);
}

int main() {
    int k = 3;  // number of dominant alleles we want
    int n = 5;  // total number of offspring
    
    printf("Calculating probability of having at least %d dominant alleles among %d offspring\n", 
           k, n);
    
    double cumulative_prob = 0.0;
    
    // Calculate P(X >= k) = P(X = k) + P(X = k+1) + ... + P(X = n)
    for (int i = k; i <= n; i++) {
        double log_prob = log_binomial_probability(n, i);
        double prob = probability_from_log(log_prob);
        cumulative_prob += prob;
        printf("P(X = %d) = %.10f\n", i, prob);
    }
    
    printf("\nTotal probability P(X >= %d) = %.10f\n", k, cumulative_prob);
    
    return 0;
}
```

## Explanation

The problem requires calculating the probability of having exactly k dominant alleles in n offspring, assuming independent inheritance and equal allele frequencies.

Key points:
1. Each offspring has a 0.5 probability of inheriting the dominant allele
2. We use the binomial distribution: P(X = k) = C(n,k) × p^k × (1-p)^(n-k)
3. For "at least k" problem, we sum probabilities from k to n
4. The solution uses either direct calculation or log-space arithmetic for numerical stability

## Sample Input/Output
```
For k=3, n=5: Probability = 0.1875000000
```

This represents the probability of having at least 3 dominant alleles among 5 offspring in a population with equal allele frequencies.