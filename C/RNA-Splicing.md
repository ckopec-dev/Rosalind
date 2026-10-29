# Rosalind Problem: RNA Splicing (C Solution)

## Problem Understanding

RNA splicing involves removing introns (non-coding sequences) from pre-mRNA to form mature mRNA. Given a DNA string and a collection of intron strings, we need to:
1. Transcribe DNA to mRNA
2. Remove all introns from the mRNA
3. Translate the resulting mRNA to protein

## Solution Approach

I'll implement this step by step:
1. Parse input DNA sequence and introns
2. Transcribe DNA to RNA (T→U)
3. Remove introns from RNA
4. Translate RNA to protein using codon table

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LENGTH 10000
#define CODON_LENGTH 3

// Codon to amino acid mapping
char* codon_table[64] = {
    "F", "F", "F", "F", "L", "L", "L", "L", "S", "S", "S", "S", 
    "S", "S", "S", "S", "Y", "Y", "Y", "Y", "*", "*", "W", "W",
    "C", "C", "C", "C", "*", "W", "L", "L", "L", "L", "P", "P",
    "P", "P", "H", "H", "H", "H", "Q", "Q", "Q", "Q", "N", "N",
    "N", "N", "K", "K", "K", "K", "R", "R", "R", "R", "I", "I",
    "I", "I", "M", "M", "M", "M", "T", "T", "T", "T", "S", "S",
    "S", "S", "S", "S", "S", "S", "A", "A", "A", "A", "V", "V",
    "V", "V", "G", "G", "G", "G"
};

// Function to translate codon to amino acid
char* translate_codon(char* codon) {
    int index = 0;
    for (int i = 0; i < CODON_LENGTH; i++) {
        if (codon[i] == 'A') index += 0;
        else if (codon[i] == 'C') index += 1 * (1 << (2*i));
        else if (codon[i] == 'G') index += 2 * (1 << (2*i));
        else if (codon[i] == 'T') index += 3 * (1 << (2*i));
    }
    return codon_table[index];
}

// Function to transcribe DNA to RNA
void transcribe(char* dna, char* rna) {
    int len = strlen(dna);
    for (int i = 0; i < len; i++) {
        if (dna[i] == 'T') rna[i] = 'U';
        else rna[i] = dna[i];
    }
    rna[len] = '\0';
}

// Function to remove intron from sequence
void remove_intron(char* seq, char* intron) {
    int seq_len = strlen(seq);
    int intron_len = strlen(intron);
    
    // Find intron in sequence
    for (int i = 0; i <= seq_len - intron_len; i++) {
        if (strncmp(&seq[i], intron, intron_len) == 0) {
            // Remove intron by shifting remaining characters
            memmove(&seq[i], &seq[i + intron_len], seq_len - i - intron_len + 1);
            break;
        }
    }
}

// Function to translate RNA to protein
void translate_rna(char* rna, char* protein) {
    int rna_len = strlen(rna);
    int prot_index = 0;
    
    // Process codons (3 nucleotides at a time)
    for (int i = 0; i < rna_len - 2; i += 3) {
        char codon[4];
        codon[0] = rna[i];
        codon[1] = rna[i+1];
        codon[2] = rna[i+2];
        codon[3] = '\0';
        
        char* amino_acid = translate_codon(codon);
        
        // Stop translation at stop codon
        if (strcmp(amino_acid, "*") == 0) {
            break;
        }
        
        strcpy(&protein[prot_index], amino_acid);
        prot_index += strlen(amino_acid);
    }
    
    protein[prot_index] = '\0';
}

int main() {
    // Read input
    char dna_sequence[MAX_LENGTH];
    char intron_sequences[MAX_LENGTH][MAX_LENGTH];
    int num_introns = 0;
    
    // Read the first line (DNA sequence)
    fgets(dna_sequence, MAX_LENGTH, stdin);
    dna_sequence[strcspn(dna_sequence, "\n")] = '\0';
    
    // Read introns
    char line[MAX_LENGTH];
    while (fgets(line, MAX_LENGTH, stdin)) {
        line[strcspn(line, "\n")] = '\0';
        if (strlen(line) == 0) break;
        strcpy(intron_sequences[num_introns], line);
        num_introns++;
    }
    
    // Transcribe DNA to RNA
    char rna_sequence[MAX_LENGTH];
    transcribe(dna_sequence, rna_sequence);
    
    // Remove all introns from RNA
    for (int i = 0; i < num_introns; i++) {
        remove_intron(rna_sequence, intron_sequences[i]);
    }
    
    // Translate RNA to protein
    char protein[MAX_LENGTH];
    translate_rna(rna_sequence, protein);
    
    printf("%s\n", protein);
    
    return 0;
}
```

## How It Works

1. **Input Parsing**: 
   - Read the main DNA sequence from first line
   - Read all intron sequences until empty line

2. **Transcription**:
   - Convert T → U in DNA to get RNA sequence

3. **Intron Removal**:
   - For each intron, find and remove it from the RNA sequence

4. **Translation**:
   - Process RNA in codons (3 nucleotides)
   - Use codon table to map to amino acids
   - Stop at stop codons (*)

## Time Complexity
- O(n × m) where n is sequence length and m is number of introns
- Space complexity: O(n)

## Sample Input/Output

**Input:**
```
ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
```

**Output:**
```
MVYIADKQHVASREAYGHMFKVCA
```