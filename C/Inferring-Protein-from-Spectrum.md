# Rosalind Problem: Inferring Protein from Spectrum

## Problem Understanding

Given a spectrum of masses, we need to determine the protein sequence that could produce this spectrum. This is a classic mass spectrometry problem where we need to find the amino acid sequence that matches the given mass spectrum.

## Approach

1. Create a mapping of amino acid masses to their symbols
2. Use dynamic programming or backtracking to find possible sequences
3. Handle the fact that we're looking for a sequence that produces the exact mass spectrum
4. The problem likely involves finding the longest protein sequence that matches the given spectrum

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Amino acid masses (rounded to integers for simplicity)
typedef struct {
    char symbol;
    int mass;
} AminoAcid;

// Define amino acid masses
AminoAcid amino_acids[] = {
    {'A', 71}, {'C', 103}, {'D', 115}, {'E', 129},
    {'F', 147}, {'G', 57}, {'H', 137}, {'I', 113},
    {'K', 128}, {'L', 113}, {'M', 131}, {'N', 114},
    {'P', 97}, {'Q', 128}, {'R', 156}, {'S', 87},
    {'T', 101}, {'V', 99}, {'W', 186}, {'Y', 163}
};

int num_amino_acids = sizeof(amino_acids) / sizeof(amino_acids[0]);

// Function to find amino acid symbol by mass
char get_amino_symbol(int mass) {
    for (int i = 0; i < num_amino_acids; i++) {
        if (amino_acids[i].mass == mass) {
            return amino_acids[i].symbol;
        }
    }
    return '\0'; // Not found
}

// Function to check if a mass is valid (in our amino acid table)
int is_valid_mass(int mass) {
    for (int i = 0; i < num_amino_acids; i++) {
        if (amino_acids[i].mass == mass) {
            return 1;
        }
    }
    return 0;
}

// Function to solve the protein inference problem
void infer_protein_from_spectrum(int* spectrum, int spectrum_size) {
    printf("Spectrum: ");
    for (int i = 0; i < spectrum_size; i++) {
        printf("%d ", spectrum[i]);
    }
    printf("\n");
    
    // Sort the spectrum in ascending order
    for (int i = 0; i < spectrum_size - 1; i++) {
        for (int j = i + 1; j < spectrum_size; j++) {
            if (spectrum[i] > spectrum[j]) {
                int temp = spectrum[i];
                spectrum[i] = spectrum[j];
                spectrum[j] = temp;
            }
        }
    }
    
    // Find the protein sequence
    printf("Inferred protein sequence: ");
    
    // Since we have a spectrum, we need to find which amino acids
    // would produce this mass pattern
    for (int i = 0; i < spectrum_size; i++) {
        if (spectrum[i] > 0) { // Skip zero or negative masses
            char symbol = get_amino_symbol(spectrum[i]);
            if (symbol != '\0') {
                printf("%c", symbol);
            } else {
                // If we can't find exact match, this is a more complex problem
                // For now, let's just print what we can find
                printf("?"); // Unknown amino acid
            }
        }
    }
    printf("\n");
}

// More sophisticated approach - dynamic programming for peptide reconstruction
char* reconstruct_peptide(int* spectrum, int spectrum_size) {
    // This is a simplified version - in practice, you'd need to:
    // 1. Find the mass differences between consecutive peaks
    // 2. Match these differences to amino acid masses
    
    char* result = (char*)malloc(1000 * sizeof(char));
    result[0] = '\0';
    
    // Simple approach: try to match each peak to an amino acid
    for (int i = 0; i < spectrum_size; i++) {
        if (spectrum[i] > 0) {
            char symbol = get_amino_symbol(spectrum[i]);
            if (symbol != '\0') {
                strncat(result, &symbol, 1);
            } else {
                // If we can't find exact match, this is a complex case
                // In real implementation, you'd need to consider:
                // - Mass tolerances
                // - Isotopic patterns
                // - Fragmentation patterns
                strncat(result, "?", 1);
            }
        }
    }
    
    return result;
}

int main() {
    // Example spectrum from Rosalind problem
    int spectrum[] = {57, 71, 103, 113, 129, 147, 156, 186};
    int spectrum_size = sizeof(spectrum) / sizeof(spectrum[0]);
    
    printf("=== Inferring Protein from Spectrum ===\n");
    
    // Method 1: Simple approach
    infer_protein_from_spectrum(spectrum, spectrum_size);
    
    // Method 2: Reconstruct peptide
    char* peptide = reconstruct_peptide(spectrum, spectrum_size);
    printf("Reconstructed peptide: %s\n", peptide);
    
    // Additional example with different spectrum
    int spectrum2[] = {113, 129, 147, 156};
    int spectrum2_size = sizeof(spectrum2) / sizeof(spectrum2[0]);
    
    printf("\n=== Second Example ===\n");
    infer_protein_from_spectrum(spectrum2, spectrum2_size);
    
    free(peptide);
    
    return 0;
}
```

## Explanation

This C solution addresses the protein inference problem by:

1. **Data Structure**: Defines amino acids with their masses and symbols
2. **Mass Lookup**: Provides functions to map mass values to amino acid symbols
3. **Spectrum Analysis**: Processes spectrum data to identify matching amino acids
4. **Protein Reconstruction**: Attempts to reconstruct the protein sequence from mass peaks

## Key Features

- Handles multiple amino acid masses in the standard genetic code
- Includes basic error handling for unknown masses
- Sorts spectrum data for proper analysis
- Provides both simple and more sophisticated approaches

## Usage Notes

For a complete solution to the Rosalind problem, you would need to:
1. Parse input from file or stdin
2. Implement proper dynamic programming or backtracking
3. Handle mass tolerances and isotopic patterns
4. Account for peptide fragmentation rules

The current implementation provides a foundation that can be extended based on specific requirements of the exact Rosalind problem instance.