# Rosalind Problem: Find Substrings of a Genome Encoding a Given Amino Acid String

## Problem Understanding

We need to find all substrings of a given DNA sequence that encode a specific amino acid sequence using the standard genetic code.

## Solution Approach

I'll implement this in Assembly language, focusing on:
1. Reading the DNA sequence and target amino acid string
2. Translating DNA codons to amino acids
3. Finding matching substrings
4. Returning valid start positions

## Assembly Implementation

```assembly
.section .data
    # Genetic code table (codon -> amino acid)
    codon_table:
        .word 0x45454545   # AAA -> K
        .word 0x45454545   # AAG -> K  
        .word 0x45454545   # AAC -> N
        .word 0x45454545   # AAT -> N
        .word 0x41414141   # AAA -> K
        .word 0x41414141   # AAG -> K
        .word 0x41414141   # AAC -> N
        .word 0x41414141   # AAT -> N
        # ... (complete genetic code table)
    
    # Example codon mappings
    .word 0x47474747   # GGG -> G
    .word 0x47474747   # GGA -> G
    .word 0x47474747   # GGC -> A
    .word 0x47474747   # GGT -> V
    
    # Target amino acid string
    target_aa: .ascii "MA"
    target_len: .long 2
    
    # DNA sequence
    dna_seq: .ascii "ATGGCCATGGCCCCCCTCGTGC"
    dna_len: .long 24

.section .text
    .global _start

_start:
    # Initialize registers
    movl dna_len(%esp), %ecx        # Length of DNA sequence
    movl $0, %esi                   # Start position counter
    movl target_len(%esp), %edi     # Length of target amino acid string
    
    # Loop through all possible substrings
check_substrings:
    cmpl %ecx, %esi                 # Check if we've gone too far
    jge done                        # If so, we're done
    
    # Calculate substring length (3 * target_aa_length)
    movl %edi, %eax
    shll $2, %eax                   # Multiply by 4 (for 4-bit encoding)
    movl %eax, %ebx                 # Store in ebx
    
    # Check if substring is long enough
    cmpl %ebx, %ecx
    jl next_position                  # Skip if not long enough
    
    # Translate current substring to amino acids
    call translate_codons
    
    # Compare with target amino acid sequence
    call compare_aa_sequences
    
    jmp next_position

next_position:
    incl %esi                       # Move to next position
    jmp check_substrings

translate_codons:
    # Input: %esi = start position in DNA
    # Output: translated amino acids in buffer
    pushl %esi
    pushl %ecx
    pushl %edi
    
    movl %esi, %eax                 # Current position
    movl $0, %ebx                   # Buffer index
    movl target_len(%esp), %ecx     # Number of codons to translate
    
translate_loop:
    cmpl $0, %ecx
    jz translate_done
    
    # Get 3 nucleotides (codon)
    movb dna_seq(%eax), %dl         # First nucleotide
    movb dna_seq+1(%eax), %dh       # Second nucleotide
    movb dna_seq+2(%eax), %cl       # Third nucleotide
    
    # Convert codon to amino acid using lookup table
    call codon_to_aa
    
    # Store result in buffer
    movb %al, translated_buffer(%ebx)
    
    addl $3, %eax                   # Move to next codon
    incl %ebx                       # Next buffer position
    decl %ecx                       # Decrement counter
    
    jmp translate_loop
    
translate_done:
    popl %edi
    popl %ecx
    popl %esi
    ret

codon_to_aa:
    # Convert 3 nucleotides to amino acid code
    # This is a simplified version - full implementation would need 
    # complete genetic code mapping
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
    
    # Simple lookup (in practice, this would be more complex)
    movb %dl, %al                   # First nucleotide
    movb %dh, %ah                   # Second nucleotide  
    movb %cl, %bl                   # Third nucleotide
    
    # Lookup in codon table
    # For demonstration purposes:
    cmpb $0x41, %dl                 # 'A' 
    je aa_A
    cmpb $0x43, %dl                 # 'C'
    je aa_C
    cmpb $0x47, %dl                 # 'G'
    je aa_G
    cmpb $0x54, %dl                 # 'T'
    je aa_T
    
aa_A:
    movb $0x4D, %al                 # 'M' (Methionine)
    jmp codon_done
    
aa_C:
    movb $0x41, %al                 # 'A' (Alanine)
    jmp codon_done
    
aa_G:
    movb $0x47, %al                 # 'G' (Glycine)
    jmp codon_done
    
aa_T:
    movb $0x56, %al                 # 'V' (Valine)
    jmp codon_done
    
codon_done:
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
    ret

compare_aa_sequences:
    # Compare translated amino acids with target
    pushl %esi
    pushl %edi
    pushl %ecx
    
    movl %esi, %eax                 # Start position
    movl target_len(%esp), %ecx     # Length of target
    
compare_loop:
    cmpl $0, %ecx
    jz match_found
    
    # Get current amino acid from translation buffer
    movb translated_buffer(%edi), %dl
    movb target_aa(%esi), %dh
    
    cmpb %dh, %dl
    jne no_match
    
    incl %esi
    incl %edi
    decl %ecx
    jmp compare_loop
    
match_found:
    # Found matching substring - record position
    pushl %eax                      # Push start position
    jmp compare_done
    
no_match:
    # No match, continue searching
    movl $0, %eax                   # No match flag
    
compare_done:
    popl %ecx
    popl %edi
    popl %esi
    ret

done:
    # Exit program
    movl $1, %eax                   # sys_exit
    movl $0, %ebx                   # exit status
    int $0x80
```

## Key Components

1. **Genetic Code Table**: Maps codons (3-nucleotide sequences) to amino acids
2. **Substring Generation**: Checks all possible substrings of appropriate length
3. **Codon Translation**: Converts DNA codons to amino acid sequences
4. **Comparison Logic**: Matches translated sequences with target amino acid string

## Time Complexity

- O(n × m) where n is the length of DNA sequence and m is the length of target amino acid string
- For each position, we translate a substring and compare it

## Space Complexity

- O(m) for storing translated amino acid buffer where m is the target length

This implementation provides the basic framework for finding substrings that encode a given amino acid sequence in Assembly language. The genetic code table would need to be fully populated with all 64 codon mappings for complete functionality.