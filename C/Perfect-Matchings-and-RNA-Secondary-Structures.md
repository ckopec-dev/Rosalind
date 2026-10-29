# Rosalind Problem: Perfect Matchings and RNA Secondary Structures

## Problem Understanding

This problem asks us to find the number of perfect matchings in an RNA secondary structure where:
- Each nucleotide can pair with its complement (A↔U, G↔C)
- A perfect matching means every nucleotide is paired exactly once
- We need to count all possible perfect matchings

## Approach

The key insight is that for a perfect matching to exist, we need an even number of each nucleotide type. If we have:
- n_A = number of A's 
- n_U = number of U's
- n_G = number of G's  
- n_C = number of C's

For a perfect matching, we need:
- n_A = n_U (A pairs with U)
- n_G = n_C (G pairs with C)

The number of perfect matchings is then the product of:
- n_A! / ((n_A/2)!)² for A-U pairings
- n_G! / ((n_G/2)!)² for G-C pairings

## Solution in C

```c
#include <stdio.h>
#include <string.h>

// Function to calculate factorial
long long factorial(int n) {
    if (n <= 1) return 1;
    long long result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

// Function to calculate combination nCk
long long combination(int n, int k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    
    // Use the smaller of k and n-k for efficiency
    if (k > n - k) k = n - k;
    
    long long result = 1;
    for (int i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

// Function to calculate number of perfect matchings
long long count_perfect_matchings(char* sequence) {
    int n = strlen(sequence);
    
    // Count nucleotides
    int a_count = 0, u_count = 0, g_count = 0, c_count = 0;
    
    for (int i = 0; i < n; i++) {
        switch (sequence[i]) {
            case 'A': a_count++; break;
            case 'U': u_count++; break;
            case 'G': g_count++; break;
            case 'C': c_count++; break;
        }
    }
    
    // For perfect matching, we need equal pairs
    // A must pair with U, G must pair with C
    if (a_count != u_count || g_count != c_count) {
        return 0; // Impossible to have perfect matching
    }
    
    // If we have n_A = n_U = a and n_G = n_C = g,
    // then the number of ways is:
    // (2*a)! / (2^a * a!) * (2*g)! / (2^g * g!)
    
    long long result = 1;
    
    // Calculate for A-U pairs
    if (a_count > 0) {
        // This is actually the double factorial approach
        // For n pairs, it's (2n)! / (2^n * n!)
        long long numerator = factorial(2 * a_count);
        long long denominator = 1;
        for (int i = 0; i < a_count; i++) {
            denominator *= 2;
        }
        denominator *= factorial(a_count);
        
        result *= (numerator / denominator);
    }
    
    // Calculate for G-C pairs
    if (g_count > 0) {
        long long numerator = factorial(2 * g_count);
        long long denominator = 1;
        for (int i = 0; i < g_count; i++) {
            denominator *= 2;
        }
        denominator *= factorial(g_count);
        
        result *= (numerator / denominator);
    }
    
    return result;
}

// Better approach using the mathematical formula
long long count_perfect_matchings_optimized(char* sequence) {
    int n = strlen(sequence);
    
    // Count nucleotides
    int a_count = 0, u_count = 0, g_count = 0, c_count = 0;
    
    for (int i = 0; i < n; i++) {
        switch (sequence[i]) {
            case 'A': a_count++; break;
            case 'U': u_count++; break;
            case 'G': g_count++; break;
            case 'C': c_count++; break;
        }
    }
    
    // For perfect matching, we need equal pairs
    if (a_count != u_count || g_count != c_count) {
        return 0; // Impossible to have perfect matching
    }
    
    // The number of perfect matchings for n pairs is (2n)!/(2^n * n!)
    // This can be computed as the product of (2i-1) for i from 1 to n
    // which equals (2n-1)!!
    
    long long result = 1;
    
    // For A-U pairs
    if (a_count > 0) {
        // This is actually the Catalan number related calculation
        // But we're looking for perfect matchings in a specific context
        // For n pairs, it's (2n-1)!! = (2n-1)*(2n-3)*...*3*1
        long long temp_result = 1;
        for (int i = 1; i <= a_count; i++) {
            temp_result *= (2 * i - 1);
        }
        result *= temp_result;
    }
    
    // For G-C pairs  
    if (g_count > 0) {
        long long temp_result = 1;
        for (int i = 1; i <= g_count; i++) {
            temp_result *= (2 * i - 1);
        }
        result *= temp_result;
    }
    
    return result;
}

// Simple and correct approach using the standard formula
long long count_perfect_matchings_correct(char* sequence) {
    int n = strlen(sequence);
    
    // Count nucleotides
    int a_count = 0, u_count = 0, g_count = 0, c_count = 0;
    
    for (int i = 0; i < n; i++) {
        switch (sequence[i]) {
            case 'A': a_count++; break;
            case 'U': u_count++; break;
            case 'G': g_count++; break;
            case 'C': c_count++; break;
        }
    }
    
    // For perfect matching, we need equal pairs
    if (a_count != u_count || g_count != c_count) {
        return 0; // Impossible to have perfect matching
    }
    
    // The number of perfect matchings in a sequence with n_A=n_U and n_G=n_C
    // is given by: (n_A)! * (n_G)! / (2^(n_A+n_G)) * (n_A+n_G)!
    // But for the specific case where we have pairs of A-U and G-C:
    // It's actually the double factorial approach
    
    long long result = 1;
    
    // For each type of pairing, calculate the number of ways
    // This is equivalent to calculating (2n-1)!! for n pairs
    if (a_count > 0) {
        long long temp = 1;
        for (int i = 1; i <= a_count; i++) {
            temp *= (2 * i - 1);
        }
        result *= temp;
    }
    
    if (g_count > 0) {
        long long temp = 1;
        for (int i = 1; i <= g_count; i++) {
            temp *= (2 * i - 1);
        }
        result *= temp;
    }
    
    return result;
}

int main() {
    char sequence[1000];
    
    printf("Enter RNA sequence: ");
    scanf("%s", sequence);
    
    long long result = count_perfect_matchings_correct(sequence);
    
    printf("Number of perfect matchings: %lld\n", result);
    
    return 0;
}
```

## Explanation

This solution works as follows:

1. **Count nucleotides**: Count how many A's, U's, G's, and C's are in the sequence
2. **Validate perfect matching**: Check if A's equal U's and G's equal C's (necessary condition)
3. **Calculate result**: For n pairs of each type, the number of perfect matchings is (2n-1)!! for each pairing type
4. **Multiply results**: Since we have independent pairings, multiply the results together

The mathematical approach uses double factorials: for n pairs, there are (2n-1)!! ways to arrange them.

## Example Usage

For sequence "AUAU":
- A count = 2, U count = 2
- Perfect matching possible
- Result = 3 × 3 = 9 (since we have 2 A-U pairs and 0 G-C pairs)

For sequence "AAUU":
- A count = 2, U count = 2  
- Perfect matching possible
- Result = 3 (for the 2 A-U pairs)