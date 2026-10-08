# Rosalind Problem: Translate an RNA String into an Amino Acid String

## Problem Understanding

We need to translate an RNA string into an amino acid string using the genetic code. Each codon (3 nucleotides) corresponds to a specific amino acid.

## Assembly Solution

```assembly
; Translate RNA string to amino acid string
; Input: RNA string in memory
; Output: Amino acid string in memory

.data
    ; Genetic code table (codon -> amino acid)
    codon_table:
        .word 'UUU', 'F'      ; Phenylalanine
        .word 'UUC', 'F'      ; Phenylalanine
        .word 'UUA', 'L'      ; Leucine
        .word 'UUG', 'L'      ; Leucine
        .word 'CUU', 'L'      ; Leucine
        .word 'CUC', 'L'      ; Leucine
        .word 'CUA', 'L'      ; Leucine
        .word 'CUG', 'L'      ; Leucine
        .word 'AUU', 'I'      ; Isoleucine
        .word 'AUC', 'I'      ; Isoleucine
        .word 'AUA', 'I'      ; Isoleucine
        .word 'UGU', 'C'      ; Cysteine
        .word 'UGC', 'C'      ; Cysteine
        .word 'UGG', 'W'      ; Tryptophan
        .word 'CCU', 'P'      ; Proline
        .word 'CCC', 'P'      ; Proline
        .word 'CCA', 'P'      ; Proline
        .word 'CCG', 'P'      ; Proline
        .word 'ACU', 'T'      ; Threonine
        .word 'ACC', 'T'      ; Threonine
        .word 'ACA', 'T'      ; Threonine
        .word 'ACG', 'T'      ; Threonine
        .word 'GCU', 'A'      ; Alanine
        .word 'GCC', 'A'      ; Alanine
        .word 'GCA', 'A'      ; Alanine
        .word 'GCG', 'A'      ; Alanine
        .word 'UAU', 'Y'      ; Tyrosine
        .word 'UAC', 'Y'      ; Tyrosine
        .word 'UGA', '*'      ; Stop codon
        .word 'UAG', '*'      ; Stop codon
        .word 'UGG', 'W'      ; Tryptophan
        .word 'CAU', 'H'      ; Histidine
        .word 'CAC', 'H'      ; Histidine
        .word 'CAA', 'Q'      ; Glutamine
        .word 'CAG', 'Q'      ; Glutamine
        .word 'AAU', 'N'      ; Asparagine
        .word 'AAC', 'N'      ; Asparagine
        .word 'AAA', 'K'      ; Lysine
        .word 'AAG', 'K'      ; Lysine
        .word 'GAU', 'D'      ; Aspartic acid
        .word 'GAC', 'D'      ; Aspartic acid
        .word 'GAA', 'E'      ; Glutamic acid
        .word 'GAG', 'E'      ; Glutamic acid
        .word 'GGU', 'G'      ; Glycine
        .word 'GGC', 'G'      ; Glycine
        .word 'GGA', 'G'      ; Glycine
        .word 'GGG', 'G'      ; Glycine
        .word 0               ; End marker

    rna_string: .ascii "AUGGCCUUU"
    rna_length: .long 9
    
    amino_acid_output: .space 100

.text
.globl _start

_start:
    ; Initialize registers
    movl rna_length, %ecx        ; Load RNA string length
    movl $rna_string, %esi       ; Load RNA string address
    movl $amino_acid_output, %edi ; Load output buffer address
    xorl %eax, %eax              ; Clear EAX for counting
    
translate_loop:
    ; Check if we've processed all codons
    cmpl $0, %ecx
    jle translate_done
    
    ; Process 3 nucleotides (1 codon)
    movb (%esi), %al             ; Load first nucleotide
    movb 1(%esi), %bl            ; Load second nucleotide
    movb 2(%esi), %cl            ; Load third nucleotide
    
    ; Create codon value for lookup
    movl %eax, %edx              ; First nucleotide to EDX
    shll $8, %edx                ; Shift left by 8 bits
    orl %ebx, %edx               ; OR second nucleotide
    shll $8, %edx                ; Shift left by 8 more bits
    orl %ecx, %edx               ; OR third nucleotide
    
    ; Look up codon in table
    call lookup_codon
    
    ; Store result
    movb %al, (%edi)
    incb %edi                    ; Move to next output position
    
    ; Advance to next codon
    addl $3, %esi                ; Skip 3 nucleotides
    subl $3, %ecx                ; Decrement codon count
    
    jmp translate_loop

translate_done:
    ; Null terminate output string
    movb $0, (%edi)
    
    ; Exit program
    movl $1, %eax                ; sys_exit
    movl $0, %ebx                ; exit status
    int $0x80

lookup_codon:
    ; Input: EDX = codon value
    ; Output: AL = amino acid character
    pushl %edx                   ; Save codon value
    
    movl $codon_table, %esi      ; Load table address
    xorl %eax, %eax              ; Clear return value
    
lookup_loop:
    movl (%esi), %ebx            ; Load codon from table
    
    cmpl $0, %ebx                ; Check for end marker
    je lookup_not_found
    
    cmpl %edx, %ebx              ; Compare with input codon
    je lookup_found
    
    addl $4, %esi                ; Move to next entry (codon + amino acid)
    jmp lookup_loop
    
lookup_found:
    movb 4(%esi), %al            ; Load amino acid character
    popl %edx                    ; Restore codon value
    ret
    
lookup_not_found:
    movb $'X', %al               ; Default to unknown amino acid
    popl %edx                    ; Restore codon value
    ret
```

## Key Features of the Solution

1. **Genetic Code Table**: Stores codons and their corresponding amino acids
2. **Codon Processing**: Reads RNA string in groups of 3 nucleotides
3. **Lookup Function**: Maps codons to amino acid characters using a hash-like approach
4. **Memory Management**: Properly handles input/output buffers
5. **Error Handling**: Includes default character for unknown codons

## Time and Space Complexity

- **Time Complexity**: O(n) where n is the length of RNA string
- **Space Complexity**: O(1) for processing (excluding input/output storage)

This solution efficiently translates RNA sequences to amino acid sequences using assembly language constructs and proper memory management.