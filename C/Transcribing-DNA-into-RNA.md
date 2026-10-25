# Rosalind Problem: Transcribing DNA into RNA (C Solution)

## Problem Description
In DNA strings, symbols 'A', 'T', 'G', and 'C' represent the nucleotides adenine, thymine, guanine, and cytosine respectively. In RNA strings, these symbols are replaced by 'A', 'U', 'G', and 'C' respectively.

Given a DNA string `s` of length at most 1000 nucleotides, return the transcribed RNA string.

## Solution

```c
#include <stdio.h>
#include <string.h>

int main() {
    char dna[1001];  // DNA string with max length 1000 + 1 for null terminator
    char rna[1001];  // RNA string
    
    // Read the DNA string
    scanf("%s", dna);
    
    int len = strlen(dna);
    
    // Transcribe DNA to RNA
    for (int i = 0; i < len; i++) {
        switch (dna[i]) {
            case 'A':
                rna[i] = 'U';
                break;
            case 'T':
                rna[i] = 'A';
                break;
            case 'G':
                rna[i] = 'C';
                break;
            case 'C':
                rna[i] = 'G';
                break;
            default:
                rna[i] = dna[i];  // Should not happen for valid DNA
        }
    }
    
    rna[len] = '\0';  // Null terminate the RNA string
    
    // Print the result
    printf("%s\n", rna);
    
    return 0;
}
```

## Explanation

1. **Input Reading**: We read the DNA string using `scanf("%s", dna)` which stores the input in the character array `dna`.

2. **String Length**: We use `strlen()` to get the length of the DNA string.

3. **Transcription Logic**: 
   - 'A' → 'U'
   - 'T' → 'A' 
   - 'G' → 'C'
   - 'C' → 'G'

4. **Output**: We print the resulting RNA string.

## Example

**Input:**
```
GATGGAACTTGACTACGTAAATT
```

**Output:**
```
GAUGGAACUUGACUACGUAAAUU
```

## Key Points

- The solution handles DNA strings of maximum length 1000 as specified
- Uses a simple switch-case statement for the transcription mapping
- Properly null-terminates the output string
- Time complexity: O(n) where n is the length of the DNA string
- Space complexity: O(n) for storing the RNA string

This solution correctly handles the DNA to RNA transcription as required by the Rosalind problem.