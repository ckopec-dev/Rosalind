# Rosalind Problem: Counting Disease Carriers

## Problem Statement
Given three positive integers representing the number of couples in a population with genotypes AA, Aa, and aa (where A is dominant and a is recessive), determine the probability that a randomly selected individual carries at least one recessive allele (i.e., has genotype Aa or aa).

## Solution Approach

```c
#include <stdio.h>

int main() {
    // Read input values
    int AA, Aa, aa;
    scanf("%d %d %d", &AA, &Aa, &aa);
    
    // Total number of couples
    int total_couples = AA + Aa + aa;
    
    // Total number of individuals (each couple produces 2 offspring)
    int total_individuals = 2 * total_couples;
    
    // Calculate probability that a randomly selected individual carries at least one recessive allele
    // This is equivalent to 1 - probability that individual is homozygous dominant (AA)
    
    // Probability of selecting an AA individual:
    // From AA couples: each produces 2 AA offspring
    // From Aa couples: each produces 1 AA offspring out of 4 possible genotypes (AA, Aa, Aa, aa)
    // From aa couples: each produces 0 AA offspring
    
    double prob_AA = 0.0;
    
    // Probability of being AA from different types of couples:
    // AA x AA: 100% AA offspring
    prob_AA += (double)AA * 2.0 / total_individuals;
    
    // Aa x Aa: 25% AA offspring (1 out of 4)
    prob_AA += (double)Aa * 1.0 / total_individuals;
    
    // aa x aa: 0% AA offspring
    // No contribution to prob_AA
    
    // Probability of carrying at least one recessive allele = 1 - P(AA)
    double result = 1.0 - prob_AA;
    
    printf("%.5f\n", result);
    
    return 0;
}
```

Wait, let me reconsider this problem more carefully based on standard genetic principles:

```c
#include <stdio.h>

int main() {
    // Read input values (number of couples with each genotype)
    int AA, Aa, aa;
    scanf("%d %d %d", &AA, &Aa, &aa);
    
    // Total number of individuals
    int total_individuals = 2 * (AA + Aa + aa);
    
    // Calculate probability that a randomly selected individual carries at least one recessive allele
    // This means either Aa or aa genotype
    
    // Expected number of each genotype in offspring:
    // AA x AA: all offspring are AA
    // AA x Aa: 50% AA, 50% Aa  
    // AA x aa: all offspring are Aa
    // Aa x Aa: 25% AA, 50% Aa, 25% aa
    // Aa x aa: 50% Aa, 50% aa
    // aa x aa: all offspring are aa
    
    double expected_Aa = 0.0;
    double expected_aa = 0.0;
    
    // From AA x AA couples (2 individuals each): all offspring are AA
    expected_Aa += 0.0;  // No Aa offspring from these
    expected_aa += 0.0;  // No aa offspring from these
    
    // From AA x Aa couples: 50% Aa offspring  
    expected_Aa += (double)Aa * 1.0;  // Each couple contributes 1 Aa individual
    expected_aa += 0.0;  // No aa offspring
    
    // From AA x aa couples: all Aa offspring
    expected_Aa += (double)aa * 1.0;  // Each couple contributes 1 Aa individual  
    expected_aa += 0.0;  // No aa offspring
    
    // From Aa x Aa couples: 50% Aa, 25% aa offspring
    expected_Aa += (double)Aa * 0.5;  // Each couple contributes 1 Aa individual
    expected_aa += (double)Aa * 0.25;  // Each couple contributes 0.5 aa individuals
    
    // From Aa x aa couples: 50% Aa, 50% aa offspring
    expected_Aa += (double)aa * 0.5;  // Each couple contributes 0.5 Aa individual
    expected_aa += (double)aa * 0.5;  // Each couple contributes 0.5 aa individual
    
    // From aa x aa couples: all aa offspring
    expected_Aa += 0.0;  // No Aa offspring
    expected_aa += (double)aa * 1.0;  // Each couple contributes 1 aa individual
    
    // Total expected number of individuals with at least one recessive allele
    double total_recessive = expected_Aa + expected_aa;
    
    // Probability = (total recessive individuals) / (total individuals)
    double probability = total_recessive / total_individuals;
    
    printf("%.5f\n", probability);
    
    return 0;
}
```

Actually, let me re-read the problem carefully and provide a cleaner solution:

```c
#include <stdio.h>

int main() {
    int AA, Aa, aa;
    scanf("%d %d %d", &AA, &Aa, &aa);
    
    // Total number of individuals in population
    int total_individuals = 2 * (AA + Aa + aa);
    
    // Count expected individuals with at least one recessive allele (Aa or aa)
    double recessive_count = 0.0;
    
    // AA x AA: all offspring are AA (0 recessive)  
    // AA x Aa: 50% Aa, 50% AA
    recessive_count += (double)Aa * 1.0;  // 1 Aa from each Aa couple
    
    // AA x aa: all offspring are Aa  
    recessive_count += (double)aa * 1.0;  // 1 Aa from each aa couple
    
    // Aa x Aa: 50% Aa, 25% aa
    recessive_count += (double)Aa * 0.5 + (double)Aa * 0.25;  // 0.5 Aa + 0.25 aa from each Aa couple
    
    // Aa x aa: 50% Aa, 50% aa
    recessive_count += (double)aa * 0.5 + (double)aa * 0.5;  // 0.5 Aa + 0.5 aa from each aa couple
    
    // aa x aa: all offspring are aa (2 aa from each aa couple)
    recessive_count += (double)aa * 1.0;  // 1 aa from each aa couple
    
    // Probability = number of recessive individuals / total individuals
    double probability = recessive_count / total_individuals;
    
    printf("%.5f\n", probability);
    
    return 0;
}
```

Let me simplify this properly one more time:

```c
#include <stdio.h>

int main() {
    int AA, Aa, aa;
    scanf("%d %d %d", &AA, &Aa, &aa);
    
    // Total number of individuals (2 per couple)
    int total_individuals = 2 * (AA + Aa + aa);
    
    // Calculate probability that a randomly selected individual carries at least one recessive allele
    double prob_recessive = 0.0;
    
    // AA x AA couples: all offspring are AA → 0 recessive individuals  
    // AA x Aa couples: 50% Aa offspring → 1 Aa per couple
    prob_recessive += (double)Aa * 1.0;  
    
    // AA x aa couples: all offspring are Aa → 1 Aa per couple
    prob_recessive += (double)aa * 1.0;
    
    // Aa x Aa couples: 50% Aa + 25% aa offspring → 0.5 Aa + 0.25 aa per couple
    prob_recessive += (double)Aa * 0.75;
    
    // Aa x aa couples: 50% Aa + 50% aa offspring → 0.5 Aa + 0.5 aa per couple  
    prob_recessive += (double)aa * 1.0;
    
    // aa x aa couples: all offspring are aa → 1 aa per couple
    prob_recessive += (double)aa * 1.0;
    
    // Probability = total recessive individuals / total individuals
    double result = prob_recessive / total_individuals;
    
    printf("%.5f\n", result);
    
    return 0;
}
```