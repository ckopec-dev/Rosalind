# Rosalind Problem: Open Reading Frames (ORF)

## Problem Description
Given a DNA string, find all open reading frames (ORFs) that could potentially code for proteins.

## Solution in C

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// DNA codon table
char codon_table[64][4] = {
    "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
    "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
    "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
    "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
    "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
    "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
    "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
    "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"
};

char amino_acid_table[64] = {
    'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S',
    'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
    'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P',
    'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R',
    'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T',
    'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R',
    'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A',
    'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'
};

// Function to translate codon to amino acid
char translate_codon(char* codon) {
    for (int i = 0; i < 64; i++) {
        if (strcmp(codon, codon_table[i]) == 0) {
            return amino_acid_table[i];
        }
    }
    return 'X'; // Unknown codon
}

// Function to reverse complement a DNA string
void reverse_complement(char* dna, char* rc) {
    int len = strlen(dna);
    for (int i = 0; i < len; i++) {
        switch (dna[len - 1 - i]) {
            case 'A': rc[i] = 'T'; break;
            case 'T': rc[i] = 'A'; break;
            case 'G': rc[i] = 'C'; break;
            case 'C': rc[i] = 'G'; break;
        }
    }
    rc[len] = '\0';
}

// Function to find ORFs in a DNA string
void find_orfs(char* dna, char** results, int* count) {
    int len = strlen(dna);
    *count = 0;
    
    // Check all 3 reading frames for forward strand
    for (int frame = 0; frame < 3; frame++) {
        for (int i = frame; i <= len - 3; i += 3) {
            // Look for start codon (ATG)
            if (dna[i] == 'A' && dna[i+1] == 'T' && dna[i+2] == 'G') {
                char protein[1000] = "";
                int j = i;
                int protein_len = 0;
                
                // Translate until stop codon
                while (j <= len - 3) {
                    char codon[4] = {dna[j], dna[j+1], dna[j+2], '\0'};
                    char amino_acid = translate_codon(codon);
                    
                    if (amino_acid == '*') break; // Stop codon
                    
                    protein[protein_len++] = amino_acid;
                    j += 3;
                }
                
                if (protein_len > 0) {
                    protein[protein_len] = '\0';
                    results[*count] = malloc((strlen(protein) + 1) * sizeof(char));
                    strcpy(results[*count], protein);
                    (*count)++;
                }
            }
        }
    }
}

int main() {
    // Read input DNA string
    char dna[10000];
    fgets(dna, sizeof(dna), stdin);
    
    // Remove newline if present
    int len = strlen(dna);
    if (dna[len-1] == '\n') {
        dna[len-1] = '\0';
        len--;
    }
    
    // Create reverse complement
    char rc[10000];
    reverse_complement(dna, rc);
    
    // Find ORFs in both strands
    char* all_orfs[1000];
    int count = 0;
    
    // Forward strand
    find_orfs(dna, all_orfs, &count);
    
    // Reverse complement strand
    find_orfs(rc, all_orfs, &count);
    
    // Print results (remove duplicates)
    char* unique_orfs[1000];
    int unique_count = 0;
    
    for (int i = 0; i < count; i++) {
        int is_duplicate = 0;
        for (int j = 0; j < unique_count; j++) {
            if (strcmp(all_orfs[i], unique_orfs[j]) == 0) {
                is_duplicate = 1;
                break;
            }
        }
        
        if (!is_duplicate) {
            unique_orfs[unique_count] = malloc((strlen(all_orfs[i]) + 1) * sizeof(char));
            strcpy(unique_orfs[unique_count], all_orfs[i]);
            unique_count++;
        }
    }
    
    // Print unique ORFs
    for (int i = 0; i < unique_count; i++) {
        printf("%s\n", unique_orfs[i]);
        free(unique_orfs[i]);
    }
    
    // Free memory
    for (int i = 0; i < count; i++) {
        free(all_orfs[i]);
    }
    
    return 0;
}
```

## Explanation

This C program solves the Open Reading Frames problem by:

1. **Reading DNA input**: Takes a DNA string from standard input
2. **Reverse complement**: Creates the reverse complement of the DNA strand
3. **Finding ORFs**: Scans both strands in all three reading frames (0, 1, 2)
4. **Translation**: Converts codons to amino acids using the genetic code
5. **Stop codon detection**: Stops translation at stop codons (*)
6. **Duplicate removal**: Eliminates duplicate protein sequences

## Key Features

- Implements the standard genetic code table
- Handles both forward and reverse complement strands
- Searches all three reading frames for each strand
- Properly manages memory allocation and deallocation
- Removes duplicate protein sequences from output

## Input Format
A single DNA string on one line

## Output Format
Each unique open reading frame (protein sequence) on a separate line

The program correctly identifies all possible ORFs that could code for proteins in both directions of the DNA sequence.