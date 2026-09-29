# Rosalind Problem: Sequence a Peptide

## Problem Understanding

This problem asks us to find all possible RNA sequences that could encode a given peptide sequence, considering the genetic code and degeneracy of codons.

## Solution Approach

1. **Genetic Code Mapping**: Create a mapping from codons to amino acids
2. **Reverse Mapping**: Create a mapping from amino acids to possible codons
3. **Dynamic Programming**: Use DP to count all valid RNA sequences
4. **Backtracking**: Generate actual sequences if needed

## Assembly Implementation

```assembly
; Rosalind Sequence a Peptide Problem
; Assembly language solution

.data
    ; Genetic code table (codon -> amino acid)
    codon_table:
        .ascii "TTT" ; F
        .byte 'F'
        .ascii "TTC" ; F
        .byte 'F'
        .ascii "TTA" ; L
        .byte 'L'
        .ascii "TTG" ; L
        .byte 'L'
        .ascii "TCT" ; S
        .byte 'S'
        .ascii "TCC" ; S
        .byte 'S'
        .ascii "TCA" ; S
        .byte 'S'
        .ascii "TCG" ; S
        .byte 'S'
        .ascii "TAT" ; Y
        .byte 'Y'
        .ascii "TAC" ; Y
        .byte 'Y'
        .ascii "TAA" ; *
        .byte '*'
        .ascii "TAG" ; *
        .byte '*'
        .ascii "TGA" ; *
        .byte '*'
        .ascii "TGT" ; C
        .byte 'C'
        .ascii "TGC" ; C
        .byte 'C'
        .ascii "TGG" ; W
        .byte 'W'
        .ascii "CTT" ; L
        .byte 'L'
        .ascii "CTC" ; L
        .byte 'L'
        .ascii "CTA" ; L
        .byte 'L'
        .ascii "CTG" ; L
        .byte 'L'
        .ascii "CCT" ; P
        .byte 'P'
        .ascii "CCC" ; P
        .byte 'P'
        .ascii "CCA" ; P
        .byte 'P'
        .ascii "CCG" ; P
        .byte 'P'
        .ascii "CAT" ; H
        .byte 'H'
        .ascii "CAC" ; H
        .byte 'H'
        .ascii "CAA" ; Q
        .byte 'Q'
        .ascii "CAG" ; Q
        .byte 'Q'
        .ascii "CGT" ; R
        .byte 'R'
        .ascii "CGC" ; R
        .byte 'R'
        .ascii "CGA" ; R
        .byte 'R'
        .ascii "CGG" ; R
        .byte 'R'
        .ascii "ATT" ; I
        .byte 'I'
        .ascii "ATC" ; I
        .byte 'I'
        .ascii "ATA" ; I
        .byte 'I'
        .ascii "ATG" ; M
        .byte 'M'
        .ascii "ACT" ; T
        .byte 'T'
        .ascii "ACC" ; T
        .byte 'T'
        .ascii "ACA" ; T
        .byte 'T'
        .ascii "ACG" ; T
        .byte 'T'
        .ascii "AAT" ; N
        .byte 'N'
        .ascii "AAC" ; N
        .byte 'N'
        .ascii "AAA" ; K
        .byte 'K'
        .ascii "AAG" ; K
        .byte 'K'
        .ascii "AGT" ; S
        .byte 'S'
        .ascii "AGC" ; S
        .byte 'S'
        .ascii "AGA" ; R
        .byte 'R'
        .ascii "AGG" ; R
        .byte 'R'
        .ascii "GTT" ; V
        .byte 'V'
        .ascii "GTC" ; V
        .byte 'V'
        .ascii "GTA" ; V
        .byte 'V'
        .ascii "GTG" ; V
        .byte 'V'
        .ascii "GCT" ; A
        .byte 'A'
        .ascii "GCC" ; A
        .byte 'A'
        .ascii "GCA" ; A
        .byte 'A'
        .ascii "GCG" ; A
        .byte 'A'
        .ascii "GAT" ; D
        .byte 'D'
        .ascii "GAC" ; D
        .byte 'D'
        .ascii "GAA" ; E
        .byte 'E'
        .ascii "GAG" ; E
        .byte 'E'
        .ascii "GGT" ; G
        .byte 'G'
        .ascii "GGC" ; G
        .byte 'G'
        .ascii "GGA" ; G
        .byte 'G'
        .ascii "GGG" ; G
        .byte 'G'

    ; Codon to amino acid mapping table (size = 64 entries)
    codon_map:
        .word 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
        .word 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
        .word 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47
        .word 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63

    ; Amino acid to codon mapping (for each amino acid)
    aa_codons:
        .byte 4   ; F: TTT, TTC
        .byte 4   ; L: TTA, TTG, CTT, CTC, CTA, CTG
        .byte 4   ; S: TCT, TCC, TCA, TCG
        .byte 2   ; Y: TAT, TAC
        .byte 3   ; H: CAT, CAC
        .byte 2   ; Q: CAA, CAG
        .byte 4   ; R: CGT, CGC, CGA, CGG
        .byte 3   ; I: ATT, ATC, ATA
        .byte 1   ; M: ATG
        .byte 4   ; T: ACT, ACC, ACA, ACG
        .byte 2   ; N: AAT, AAC
        .byte 2   ; K: AAA, AAG
        .byte 4   ; V: GTT, GTC, GTA, GTG
        .byte 4   ; A: GCT, GCC, GCA, GCG
        .byte 2   ; D: GAT, GAC
        .byte 2   ; E: GAA, GAG
        .byte 4   ; G: GGT, GGC, GGA, GGG
        .byte 1   ; W: TGG
        .byte 2   ; C: TGT, TGC
        .byte 1   ; *: TAA, TAG, TGA

    peptide_input: .ascii "MA"     ; Example peptide sequence
    peptide_len: .byte 2           ; Length of peptide

    result_count: .word 0          ; To store final count

.code
main:
    ; Initialize data structures
    mov eax, offset peptide_input
    mov ecx, peptide_len
    call count_rna_sequences
    
    ; Print result
    mov eax, result_count
    call print_number
    
    ret

; Function to count all possible RNA sequences for given peptide
count_rna_sequences:
    push ebp
    mov ebp, esp
    
    ; Load parameters
    mov esi, eax        ; peptide string
    mov ecx, edx        ; peptide length
    
    ; Initialize DP array
    xor edi, edi        ; i = 0
    mov ebx, 1          ; dp[0] = 1 (empty sequence)
    
count_loop:
    cmp edi, ecx
    jge count_done
    
    ; Get current amino acid
    mov al, [esi + edi]
    
    ; Find number of codons for this amino acid
    call get_codon_count_for_aa
    
    ; Multiply to DP array
    imul ebx, eax       ; multiply with current count
    
    ; Store in DP array (simplified)
    mov [ebp-4], ebx    ; store result
    
    inc edi
    jmp count_loop

count_done:
    mov eax, ebx        ; return count
    pop ebp
    ret

; Get number of codons for given amino acid
get_codon_count_for_aa:
    push ebp
    mov ebp, esp
    
    mov al, [esi]       ; get amino acid
    sub al, 'A'         ; convert to index (assuming A=0)
    
    ; Lookup in aa_codons table
    mov bl, [aa_codons + eax]
    
    mov eax, ebx        ; return count
    pop ebp
    ret

; Helper function to print number
print_number:
    push ebp
    mov ebp, esp
    
    ; Simple conversion and printing logic
    ; In a real implementation, this would convert number to string
    ; and output it
    
    pop ebp
    ret

end main
```

## Explanation

This assembly solution:

1. **Data Section**: 
   - Defines genetic code mappings for codons to amino acids
   - Contains lookup tables for codon degeneracy
   - Stores input peptide sequence

2. **Main Algorithm**:
   - Uses dynamic programming approach where `dp[i]` represents number of ways to encode first `i` amino acids
   - For each amino acid, multiplies by the number of possible codons that can encode it
   - Returns final count

3. **Key Functions**:
   - `count_rna_sequences`: Main DP function
   - `get_codon_count_for_aa`: Looks up how many codons correspond to each amino acid

## Time Complexity
O(n) where n is the length of peptide sequence

## Space Complexity  
O(1) additional space (excluding input/output)

Note: This is a simplified implementation. A full solution would need more complex string handling and proper assembly syntax for the specific assembler being used.