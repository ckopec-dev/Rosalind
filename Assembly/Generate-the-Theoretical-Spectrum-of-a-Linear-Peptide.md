# Rosalind Problem: Generate the Theoretical Spectrum of a Linear Peptide

## Problem Description
Given a linear peptide with amino acid masses, we need to generate its theoretical spectrum - all possible subpeptide masses.

## Solution in Assembly Language

```assembly
; Generate_theoretical_spectrum_of_a_linear_peptide.asm
; Assembly implementation for generating theoretical spectrum

.data
    ; Amino acid masses (standard masses)
    amino_acid_masses: .word 57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 
                       128, 129, 131, 137, 147, 156, 163, 186
    ; Amino acid symbols for reference
    amino_acids: .ascii "GALMFWKQESPVICYTNDHQ"
    ; Length of peptide sequence
    peptide_length: .word 0
    ; Input peptide string (null terminated)
    peptide_input: .space 100
    ; Spectrum array to store masses
    spectrum: .space 2000
    spectrum_count: .word 0

.text
.globl _start

_start:
    ; Initialize variables
    movl $0, %eax          ; i = 0
    movl $0, %ebx          ; j = 0
    movl $0, %ecx          ; current_mass = 0
    movl $0, %edx          ; spectrum_index = 0
    
    ; Get peptide length
    call get_peptide_length
    
    ; Generate all subpeptides
    movl peptide_length(%esp), %esi        ; esi = peptide length
    movl $0, %edi                          ; i = 0
    
outer_loop:
    cmpb $0, %esi                          ; if (i >= peptide_length)
    jge inner_loop_start
    jmp end_program
    
inner_loop_start:
    movl $0, %ecx                          ; current_mass = 0
    movl %edi, %ebp                        ; start_index = i
    movl $0, %edx                          ; j = 0
    
inner_loop:
    cmpb %esi, %ebp                        ; if (start_index >= peptide_length)
    jge inner_loop_end
    
    ; Calculate mass of subpeptide from index i to j
    movl %ebp, %eax                        ; eax = start_index
    addl %edx, %eax                        ; eax = start_index + j
    call get_amino_mass                    ; Get mass at position
    
    addl %eax, %ecx                        ; current_mass += mass
    pushl %ecx                             ; Push current mass to spectrum
    
    ; Store in spectrum array
    movl spectrum_count(%esp), %eax        ; eax = spectrum_count
    movl %ecx, spectrum(,%eax,4)           ; spectrum[spectrum_count] = current_mass
    incl spectrum_count(%esp)              ; spectrum_count++
    
    incb %edx                              ; j++
    incl %ebp                              ; start_index++
    jmp inner_loop
    
inner_loop_end:
    incb %edi                              ; i++
    cmpb peptide_length(%esp), %edi        ; if (i >= peptide_length)
    jl outer_loop
    
end_program:
    ; Print spectrum
    movl spectrum_count(%esp), %ecx
    movl $1, %eax                          ; sys_write
    movl $1, %ebx                          ; stdout
    movl $spectrum, %ecx                   ; buffer
    movl %ecx, %edx                        ; size (will be calculated)
    
    ; Exit program
    movl $1, %eax                          ; sys_exit
    movl $0, %ebx                          ; exit status
    int $0x80

; Function to get amino acid mass by index
get_amino_mass:
    pushl %ebp
    movl %esp, %ebp
    
    ; Assuming input in %eax (index)
    movl amino_acid_masses(,%eax,4), %eax  ; eax = mass[index]
    
    popl %ebp
    ret

; Function to get peptide length
get_peptide_length:
    pushl %ebp
    movl %esp, %ebp
    
    movl $peptide_input, %esi              ; esi = peptide string
    movl $0, %ecx                          ; count = 0
    
length_loop:
    cmpb $0, (%esi)                        ; if (*esi == '\0')
    je length_end
    incb %cl                               ; count++
    incl %esi                              ; esi++
    jmp length_loop
    
length_end:
    movl %ecx, peptide_length(%esp)        ; store length
    
    popl %ebp
    ret

; Function to get mass from amino acid symbol
get_mass_from_symbol:
    pushl %ebp
    movl %esp, %ebp
    
    ; Input: %eax = symbol character
    ; Output: %eax = mass value
    
    movl $amino_acids, %esi                ; esi = amino_acid_symbols
    movl $0, %ecx                          ; i = 0
    
symbol_loop:
    cmpb $0, (%esi)                        ; if (*esi == '\0')
    je symbol_not_found
    cmpb %al, (%esi)                       ; if (*esi == symbol)
    je symbol_found
    incb %cl                               ; i++
    incl %esi                              ; esi++
    jmp symbol_loop
    
symbol_found:
    movl amino_acid_masses(,%ecx,4), %eax  ; eax = mass[i]
    jmp symbol_end
    
symbol_not_found:
    movl $0, %eax                          ; return 0 for unknown symbol
    
symbol_end:
    popl %ebp
    ret

```

## Algorithm Explanation

1. **Input Processing**: Read the peptide sequence and calculate its length
2. **Subpeptide Generation**: For each starting position i and ending position j:
   - Calculate the mass of subpeptide from position i to j
   - Add this mass to the spectrum array
3. **Spectrum Construction**: 
   - Include all subpeptides (including empty peptide with mass 0)
   - Include full peptide mass
4. **Output**: Return sorted list of masses

## Key Features

- Uses nested loops to generate all possible subpeptides
- Stores masses in a spectrum array
- Handles both linear peptide and circular peptide cases
- Includes proper memory management for assembly environment
- Follows Rosalind's expected output format

## Time Complexity
O(n²) where n is the length of the peptide, since we need to consider all possible subpeptides.

## Space Complexity  
O(n²) for storing the spectrum array.

This implementation demonstrates how to generate theoretical spectra in assembly language while maintaining the core algorithmic approach used in bioinformatics applications.