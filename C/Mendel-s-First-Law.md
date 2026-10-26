# Rosalind Problem: Mendel's First Law (MEND)

## Problem Statement
Given three positive integers k, m, and n representing the number of organisms of each genotype in a population:
- k: homozygous dominant individuals (AA)
- m: heterozygous individuals (Aa) 
- n: homozygous recessive individuals (aa)

Return the probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (i.e., the probability of getting an offspring with genotype AA, Aa, or Aa).

## Solution Approach
We need to calculate the probability that a randomly chosen pair of organisms will produce offspring with at least one dominant allele (A).

The key is to consider all possible mating pairs and their probabilities:

1. **AA × AA**: Always produces AA (probability = 1)
2. **AA × Aa**: Always produces A_ (probability = 1) 
3. **AA × aa**: Always produces Aa (probability = 1)
4. **Aa × Aa**: Produces AA (0.25), Aa (0.5), aa (0.25)
5. **Aa × aa**: Produces Aa (0.5), aa (0.5)
6. **aa × aa**: Always produces aa (probability = 0)

## C Implementation

```c
#include <stdio.h>

int main() {
    int k, m, n;
    
    // Read input values
    printf("Enter k (homozygous dominant): ");
    scanf("%d", &k);
    printf("Enter m (heterozygous): ");
    scanf("%d", &m);
    printf("Enter n (homozygous recessive): ");
    scanf("%d", &n);
    
    // Total population
    int total = k + m + n;
    
    // Calculate probability of each mating pair
    double prob = 0.0;
    
    // AA × AA
    prob += ((double)k / total) * ((double)(k - 1) / (total - 1));
    
    // AA × Aa
    prob += ((double)k / total) * ((double)m / (total - 1));
    prob += ((double)m / total) * ((double)k / (total - 1));
    
    // AA × aa
    prob += ((double)k / total) * ((double)n / (total - 1));
    prob += ((double)n / total) * ((double)k / (total - 1));
    
    // Aa × Aa
    prob += ((double)m / total) * ((double)(m - 1) / (total - 1)) * 0.75;
    
    // Aa × aa
    prob += ((double)m / total) * ((double)n / (total - 1)) * 0.5;
    prob += ((double)n / total) * ((double)m / (total - 1)) * 0.5;
    
    // aa × aa
    // This produces only recessive offspring, so it contributes 0 to our probability
    
    printf("Probability of dominant allele in offspring: %.5f\n", prob);
    
    return 0;
}
```

## Alternative Cleaner Implementation

```c
#include <stdio.h>

int main() {
    int k, m, n;
    scanf("%d %d %d", &k, &m, &n);
    
    int total = k + m + n;
    
    // Calculate all possible matings that result in dominant phenotype
    double p_dominant = 0.0;
    
    // All AA × (any) combinations - always produce dominant
    p_dominant += (double)(k * (k - 1)) / (total * (total - 1)); // AA × AA
    p_dominant += 2.0 * (double)(k * m) / (total * (total - 1));   // AA × Aa
    
    // All Aa × Aa combinations that produce dominant offspring (75%)
    p_dominant += (double)(m * (m - 1)) / (total * (total - 1)) * 0.75;
    
    // All Aa × aa combinations that produce dominant offspring (50%)
    p_dominant += 2.0 * (double)(m * n) / (total * (total - 1)) * 0.5;
    
    printf("%.5f\n", p_dominant);
    
    return 0;
}
```

## Example
Input: `k=2, m=2, n=2`
Output: `0.78333`

The solution accounts for all possible mating pairs and their probabilities of producing offspring with dominant alleles.