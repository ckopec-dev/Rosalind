# Rosalind Problem: Sex-Linked_Inheritance

## Problem Description
In sex-linked inheritance, genes located on the X chromosome are inherited differently in males and females. In this problem, we need to calculate the probability that a trait will be expressed in a population given certain conditions.

## Solution Approach
Let's break down the problem:
1. We have a population with known sex ratios
2. We need to calculate probabilities of trait expression
3. This involves understanding X-linked inheritance patterns

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

// Function to calculate probability of trait expression
double calculate_sex_linked_probability(double male_freq, double female_freq) {
    // For X-linked inheritance:
    // - Males have one X chromosome (from mother)
    // - Females have two X chromosomes (one from each parent)
    
    // If we're calculating probability that a randomly selected individual expresses the trait
    // This would depend on the specific genetic model
    
    // Simple model: assume trait is recessive on X chromosome
    // Let p be frequency of recessive allele in population
    
    // For males (XY): probability of expressing trait = p
    // For females (XX): probability of expressing trait = p^2
    
    // Overall probability = (male_freq * p) + (female_freq * p^2)
    
    // Since we don't have specific p values, let's implement a general approach
    return 0.0; // Placeholder - need more specific problem details
}

// More comprehensive implementation for typical sex-linked inheritance problems
double calculate_x_linked_probabilities(int total_males, int total_females, 
                                       double recessive_freq) {
    // In X-linked recessive inheritance:
    // - Males (XY): express trait if they inherit recessive allele from mother
    // - Females (XX): express trait only if they inherit recessive allele from both parents
    
    double male_prob = recessive_freq;  // Males have 1 X chromosome
    double female_prob = recessive_freq * recessive_freq;  // Females need 2 recessive alleles
    
    double total_prob = (total_males * male_prob + total_females * female_prob) / 
                       (total_males + total_females);
    
    return total_prob;
}

int main() {
    // Example input: number of males and females, frequency of recessive allele
    int num_males = 1000;
    int num_females = 1000;
    double freq_recessive = 0.1;  // 10% frequency of recessive allele
    
    printf("Sex-Linked Inheritance Problem Solution\n");
    printf("======================================\n");
    
    double probability = calculate_x_linked_probabilities(num_males, num_females, 
                                                        freq_recessive);
    
    printf("Number of males: %d\n", num_males);
    printf("Number of females: %d\n", num_females);
    printf("Frequency of recessive allele: %.3f\n", freq_recessive);
    printf("Probability of expressing trait in population: %.4f\n", probability);
    
    // Additional calculations
    double male_expr_prob = freq_recessive;
    double female_expr_prob = freq_recessive * freq_recessive;
    
    printf("\nDetailed probabilities:\n");
    printf("- Male probability of expressing trait: %.4f\n", male_expr_prob);
    printf("- Female probability of expressing trait: %.4f\n", female_expr_prob);
    
    return 0;
}
```

## Alternative Solution for Specific Problem Format

If the problem involves specific input/output format like in Rosalind, here's a more targeted approach:

```c
#include <stdio.h>

int main() {
    // Read input values (example format)
    int males, females;
    double p_recessive;
    
    // Example: reading from stdin or hardcoding for demonstration
    males = 500;
    females = 500;
    p_recessive = 0.25;
    
    // Calculate probabilities for X-linked inheritance
    double male_trait_prob = p_recessive;           // Males express if they have recessive allele
    double female_trait_prob = p_recessive * p_recessive;  // Females need two recessive alleles
    
    // Overall probability in population
    double total_prob = (males * male_trait_prob + females * female_trait_prob) / 
                       (males + females);
    
    printf("%.4f\n", total_prob);
    
    return 0;
}
```

## Key Concepts

1. **X-linked Inheritance**: Genes on X chromosome
2. **Males (XY)**: Have one X chromosome, so if they inherit recessive allele, they express trait
3. **Females (XX)**: Need two copies of recessive allele to express trait
4. **Population Probability**: Weighted average based on sex ratios

## Input/Output Format

For Rosalind problems, typical input would be:
```
Input: 
- Number of males and females in population
- Frequency of recessive allele (or other genetic parameters)

Output:
- Probability that an individual expresses the trait
```

The solution handles the mathematical principles of X-linked inheritance while providing a framework that can be adapted for specific problem requirements.