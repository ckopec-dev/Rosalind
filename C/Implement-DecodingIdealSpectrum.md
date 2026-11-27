# Rosalind Problem: Implement Decoding Ideal Spectrum

## Problem Understanding

This problem asks us to decode an ideal spectrum by finding the amino acid sequence that would produce the given spectrum. In an ideal spectrum, we have the masses of all possible subpeptides (prefixes and suffixes) of a peptide.

## Approach

1. We need to find a sequence of amino acids whose prefix sums match the given spectrum
2. This is essentially finding a path in a graph where edges represent amino acid masses
3. We can use dynamic programming or backtracking approach
4. The key insight is that we're looking for a sequence where cumulative sums match the spectrum values

## Solution Implementation

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Amino acid masses
const double amino_acid_masses[] = {
    57.0, 71.0, 87.0, 97.0, 99.0, 101.0, 103.0, 113.0, 114.0, 115.0,
    128.0, 129.0, 131.0, 137.0, 147.0, 156.0, 163.0, 186.0
};

const char amino_acid_codes[] = {
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V'
};

const int num_amino_acids = 18;

// Function to check if a value is approximately equal to a mass
int is_equal(double a, double b, double epsilon) {
    return (a - b) < epsilon && (b - a) < epsilon;
}

// Function to find the amino acid index for a given mass
int find_amino_acid_index(double mass, double epsilon) {
    for (int i = 0; i < num_amino_acids; i++) {
        if (is_equal(mass, amino_acid_masses[i], epsilon)) {
            return i;
        }
    }
    return -1;
}

// Function to decode the ideal spectrum
char* decode_spectrum(double* spectrum, int spectrum_size) {
    // Sort spectrum in ascending order
    for (int i = 0; i < spectrum_size - 1; i++) {
        for (int j = i + 1; j < spectrum_size; j++) {
            if (spectrum[i] > spectrum[j]) {
                double temp = spectrum[i];
                spectrum[i] = spectrum[j];
                spectrum[j] = temp;
            }
        }
    }
    
    // Create a result string
    char* result = (char*)malloc(spectrum_size * sizeof(char));
    int result_index = 0;
    
    // First element should be 0 for the empty prefix
    double current_mass = 0.0;
    double epsilon = 0.01; // Tolerance for floating point comparison
    
    // For each mass in spectrum, find which amino acid contributes to it
    for (int i = 1; i < spectrum_size; i++) {
        double diff = spectrum[i] - current_mass;
        
        int index = find_amino_acid_index(diff, epsilon);
        if (index != -1) {
            result[result_index++] = amino_acid_codes[index];
            current_mass = spectrum[i];
        } else {
            // If we can't find a matching amino acid mass, this indicates an error
            printf("Error: No matching amino acid found for mass %f\n", diff);
            free(result);
            return NULL;
        }
    }
    
    result[result_index] = '\0';
    return result;
}

// Main function to solve the problem
int main() {
    // Example spectrum from Rosalind problem
    double spectrum[] = {0.0, 113.0, 128.0, 147.0, 200.0, 215.0, 234.0, 249.0, 268.0, 283.0, 302.0, 317.0};
    int spectrum_size = 12;
    
    printf("Input spectrum: ");
    for (int i = 0; i < spectrum_size; i++) {
        printf("%.1f ", spectrum[i]);
    }
    printf("\n");
    
    char* peptide = decode_spectrum(spectrum, spectrum_size);
    
    if (peptide != NULL) {
        printf("Decoded peptide: %s\n", peptide);
        free(peptide);
    } else {
        printf("Failed to decode the spectrum.\n");
    }
    
    return 0;
}
```

## Explanation

This solution works by:

1. **Sorting** the spectrum values in ascending order
2. **Iterating** through the sorted spectrum starting from the second element (since first is always 0)
3. **Calculating differences** between consecutive elements to get the mass increments
4. **Matching** these differences with known amino acid masses
5. **Building** the peptide sequence character by character

## Key Points

- The first value in the spectrum should be 0 (representing empty prefix)
- Each subsequent value represents the cumulative mass of the peptide up to that point
- We use a tolerance of 0.01 for floating-point comparisons due to potential precision issues
- The solution assumes that we're working with the standard amino acid masses

## Time Complexity
O(n²) where n is the size of the spectrum, due to sorting and nested loops for matching masses.

## Space Complexity
O(n) for storing the result peptide sequence.

Note: This implementation handles a specific example. For the exact Rosalind problem, you would need to read input from file or standard input according to the problem requirements.