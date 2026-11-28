# Rosalind Problem: Implement PSMSearch in C

## Problem Description
Given a collection of spectra and a collection of peptides, find all peptides that could produce the given spectra through mass spectrometry.

## Solution Approach
I'll implement a peptide spectrum matching algorithm that:
1. Computes theoretical spectra for each peptide
2. Compares these with experimental spectra
3. Returns matching peptides

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PEPTIDE_LENGTH 50
#define MAX_SPECTRA_SIZE 1000
#define TOLERANCE 0.01
#define AMINO_ACID_MASS_TABLE_SIZE 26

// Mass of amino acids (rounded to 3 decimal places)
double amino_acid_mass[AMINO_ACID_MASS_TABLE_SIZE] = {
    71.037,  // A
    103.009, // C
    115.027, // D
    129.043, // E
    147.035, // F
    57.021,  // G
    137.058, // H
    113.084, // I
    128.095, // K
    138.062, // L
    113.084, // M
    131.040, // N
    87.032,  // P
    156.101, // Q
    163.063, // R
    86.037,  // S
    155.067, // T
    186.079, // V
    101.048, // W
    119.058, // Y
    0,       // Z (not used)
    0,       // X (not used)
    0,       // B (not used)
    0,       // J (not used)
    0,       // U (not used)
    0        // O (not used)
};

// Function to get mass of an amino acid
double get_amino_acid_mass(char amino_acid) {
    if (amino_acid >= 'A' && amino_acid <= 'Z') {
        return amino_acid_mass[amino_acid - 'A'];
    }
    return 0.0;
}

// Function to compute theoretical spectrum of a peptide
void compute_theoretical_spectrum(char* peptide, double* spectrum, int* spectrum_size) {
    int len = strlen(peptide);
    *spectrum_size = 0;
    
    // Add mass of water (18.0105 Da) for N-terminal
    double total_mass = 18.0105;
    
    // Compute prefix masses (N-terminal to C-terminal)
    for (int i = 0; i < len; i++) {
        total_mass += get_amino_acid_mass(peptide[i]);
        spectrum[(*spectrum_size)++] = total_mass;
    }
    
    // Add mass of water (18.0105 Da) for C-terminal
    total_mass = 18.0105;
    
    // Compute suffix masses (C-terminal to N-terminal)
    for (int i = len - 1; i >= 0; i--) {
        total_mass += get_amino_acid_mass(peptide[i]);
        spectrum[(*spectrum_size)++] = total_mass;
    }
}

// Function to check if two masses match within tolerance
int masses_match(double mass1, double mass2) {
    return fabs(mass1 - mass2) <= TOLERANCE;
}

// Function to compare theoretical spectrum with experimental spectrum
int is_spectrum_match(double* theoretical_spectrum, int theo_size, 
                      double* experimental_spectrum, int exp_size) {
    // Sort theoretical spectrum (for efficient matching)
    for (int i = 0; i < theo_size - 1; i++) {
        for (int j = i + 1; j < theo_size; j++) {
            if (theoretical_spectrum[i] > theoretical_spectrum[j]) {
                double temp = theoretical_spectrum[i];
                theoretical_spectrum[i] = theoretical_spectrum[j];
                theoretical_spectrum[j] = temp;
            }
        }
    }
    
    // Sort experimental spectrum
    for (int i = 0; i < exp_size - 1; i++) {
        for (int j = i + 1; j < exp_size; j++) {
            if (experimental_spectrum[i] > experimental_spectrum[j]) {
                double temp = experimental_spectrum[i];
                experimental_spectrum[i] = experimental_spectrum[j];
                experimental_spectrum[j] = temp;
            }
        }
    }
    
    // Match each theoretical peak with experimental peaks
    int matches = 0;
    for (int i = 0; i < theo_size; i++) {
        for (int j = 0; j < exp_size; j++) {
            if (masses_match(theoretical_spectrum[i], experimental_spectrum[j])) {
                matches++;
                break;
            }
        }
    }
    
    // Return true if at least 80% of theoretical peaks match
    return (double)matches / (double)theo_size >= 0.8;
}

// Main function to solve PSM search problem
int main() {
    // Read peptides
    char peptides[100][MAX_PEPTIDE_LENGTH];
    int peptide_count = 0;
    
    printf("Enter number of peptides: ");
    scanf("%d", &peptide_count);
    getchar(); // consume newline
    
    for (int i = 0; i < peptide_count; i++) {
        fgets(peptides[i], MAX_PEPTIDE_LENGTH, stdin);
        peptides[i][strcspn(peptides[i], "\n")] = 0; // Remove newline
    }
    
    // Read number of spectra
    int spectrum_count;
    printf("Enter number of spectra: ");
    scanf("%d", &spectrum_count);
    
    // Process each spectrum
    for (int s = 0; s < spectrum_count; s++) {
        // Read experimental spectrum (as space-separated masses)
        double experimental_spectrum[MAX_SPECTRA_SIZE];
        int exp_size = 0;
        
        printf("Enter experimental spectrum %d (space-separated values): ", s + 1);
        char line[1000];
        fgets(line, sizeof(line), stdin);
        
        char* token = strtok(line, " \n");
        while (token != NULL) {
            experimental_spectrum[exp_size++] = atof(token);
            token = strtok(NULL, " \n");
        }
        
        // For each peptide, compute theoretical spectrum and check match
        for (int i = 0; i < peptide_count; i++) {
            double theoretical_spectrum[MAX_SPECTRA_SIZE];
            int theo_size;
            
            compute_theoretical_spectrum(peptides[i], theoretical_spectrum, &theo_size);
            
            if (is_spectrum_match(theoretical_spectrum, theo_size, 
                                 experimental_spectrum, exp_size)) {
                printf("Peptide %s matches spectrum %d\n", peptides[i], s + 1);
            }
        }
    }
    
    return 0;
}
```

## Key Features of the Implementation

1. **Amino Acid Mass Table**: Stores masses of standard amino acids
2. **Theoretical Spectrum Generation**: Computes prefix and suffix masses for each peptide
3. **Spectrum Matching**: Compares theoretical and experimental spectra with tolerance
4. **80% Match Threshold**: Requires at least 80% of theoretical peaks to match

## How to Use

1. Compile: `gcc -o psm_search psm_search.c -lm`
2. Run: `./psm_search`
3. Input format:
   - Number of peptides
   - List of peptides (one per line)
   - Number of spectra
   - For each spectrum, input space-separated mass values

## Time Complexity
- O(n × m × k) where n is number of peptides, m is average peptide length, k is average spectrum size
- The matching process uses a simple nested loop approach for demonstration

## Space Complexity
- O(m + k) for storing theoretical and experimental spectra

This implementation provides a basic but functional solution to the Peptide Spectrum Matching problem as required by Rosalind's PSM Search challenge.