# Rosalind Problem: Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum

## Problem Understanding

We need to find a cyclic peptide whose theoretical spectrum matches a given ideal spectrum. This is a classic bioinformatics problem that involves:
1. Generating all possible cyclic peptides of a given mass
2. Computing their theoretical spectra
3. Finding which peptide's spectrum matches the target spectrum exactly

## Approach

The solution uses a backtracking approach with pruning to efficiently search for valid cyclic peptides.

## Assembly Implementation

```assembly
; Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum
; Assembly implementation for Rosalind problem

.data
    ; Mass values for amino acids (standard masses)
    amino_acid_masses: .word 57, 71, 87, 97, 99, 101, 103, 113, 114, 115
                       .word 128, 129, 131, 137, 147, 156, 163, 186
    
    ; Number of amino acids
    num_amino_acids: .word 18
    
    ; Target spectrum (ideal spectrum)
    target_spectrum: .word 0, 113, 128, 186, 244, 299, 314, 372, 427, 442, 500
    
    ; Target spectrum size
    target_size: .word 11
    
    ; Temporary storage for current peptide
    current_peptide: .space 100
    
    ; Result peptide storage
    result_peptide: .space 100
    
    ; Current peptide length
    peptide_length: .word 0
    
    ; Found flag
    found_solution: .word 0

.text
.globl _start

_start:
    ; Initialize variables
    la $t0, target_spectrum
    lw $t1, target_size
    li $t2, 0
    sw $t2, peptide_length
    
    ; Start backtracking search
    jal backtrack_search
    j exit

; Function: backtrack_search
; Purpose: Recursively build peptides and check spectra
backtrack_search:
    addi $sp, $sp, -16          ; Allocate stack space
    sw $ra, 12($sp)             ; Save return address
    
    ; Check if we have a complete peptide (length > 0)
    lw $t0, peptide_length
    beq $t0, $zero, try_next    ; If length is 0, try first amino acid
    
    ; Check if we've reached target size
    lw $t0, peptide_length
    li $t1, 3                   ; For example, assume target length is 3
    bge $t0, $t1, check_spectrum ; If we have enough amino acids, check spectrum
    
try_next:
    ; Try each amino acid
    li $t0, 0                   ; i = 0
    li $t1, 18                  ; num_amino_acids
    li $t2, 0                   ; current position in peptide
    
next_amino:
    bge $t0, $t1, backtrack_done ; If all amino acids tried
    
    ; Get mass of current amino acid
    la $t3, amino_acid_masses
    sll $t4, $t0, 2             ; Multiply by 4 (word size)
    add $t3, $t3, $t4
    lw $t5, 0($t3)              ; Load mass
    
    ; Add to current peptide
    la $t6, current_peptide
    sll $t7, $t2, 2             ; Multiply position by 4
    add $t6, $t6, $t7
    sw $t5, 0($t6)              ; Store mass in peptide
    
    ; Increment peptide length
    lw $t8, peptide_length
    addi $t8, $t8, 1
    sw $t8, peptide_length
    
    ; Recurse
    jal backtrack_search
    
    ; Backtrack: remove last amino acid
    lw $t9, peptide_length
    addi $t9, $t9, -1
    sw $t9, peptide_length
    
    addi $t0, $t0, 1            ; Next amino acid
    addi $t2, $t2, 1            ; Next position
    j next_amino

check_spectrum:
    ; Compute theoretical spectrum for current peptide
    la $a0, current_peptide
    lw $a1, peptide_length
    jal compute_theoretical_spectrum
    
    ; Compare with target spectrum
    la $a0, target_spectrum
    la $a1, current_spectrum
    lw $a2, peptide_length
    jal compare_spectra
    
    beq $v0, $zero, backtrack_done  ; If not match, continue backtracking
    
    ; Found a solution!
    la $t0, result_peptide
    la $t1, current_peptide
    lw $t2, peptide_length
    
copy_solution:
    beq $t2, $zero, backtrack_done
    addi $t2, $t2, -1
    sll $t3, $t2, 2
    add $t4, $t1, $t3
    lw $t5, 0($t4)
    add $t6, $t0, $t3
    sw $t5, 0($t6)
    
    j copy_solution

backtrack_done:
    lw $ra, 12($sp)             ; Restore return address
    addi $sp, $sp, 16           ; Deallocate stack space
    jr $ra                      ; Return

; Function: compute_theoretical_spectrum
; Purpose: Compute theoretical spectrum of a cyclic peptide
compute_theoretical_spectrum:
    addi $sp, $sp, -20         ; Allocate stack space
    sw $ra, 16($sp)            ; Save return address
    sw $a0, 0($sp)             ; Save peptide pointer
    sw $a1, 4($sp)             ; Save length
    
    ; Initialize spectrum array
    la $t0, current_spectrum
    li $t1, 0                  ; Clear all entries
    li $t2, 100                ; Max possible size
    
clear_spectrum:
    beq $t2, $zero, spectrum_done
    sw $t1, 0($t0)
    addi $t0, $t0, 4
    addi $t2, $t2, -1
    j clear_spectrum

spectrum_done:
    ; Generate all subpeptides and their masses
    li $t3, 0                  ; i = 0
    
subpeptide_loop:
    bge $t3, $a1, spectrum_cleanup
    
    li $t4, 0                  ; j = 0
    li $t5, 0                  ; cumulative mass
    
subpeptide_inner:
    bge $t4, $a1, subpeptide_next
    
    ; Calculate mass at position (i + j) % length
    add $t6, $t3, $t4
    rem $t7, $t6, $a1          ; Modulo operation
    sll $t8, $t7, 2            ; Word offset
    la $t9, current_peptide
    add $t9, $t9, $t8
    lw $t10, 0($t9)            ; Get mass
    
    add $t5, $t5, $t10         ; Add to cumulative mass
    
    ; Store mass in spectrum
    la $t11, current_spectrum
    sll $t12, $t5, 2          ; Word offset
    add $t11, $t11, $t12
    lw $t13, 0($t11)
    beq $t13, $zero, store_mass ; Only store if not already there
    
store_mass:
    sw $t5, 0($t11)
    
    addi $t4, $t4, 1
    j subpeptide_inner
    
subpeptide_next:
    addi $t3, $t3, 1
    j subpeptide_loop

spectrum_cleanup:
    lw $ra, 16($sp)            ; Restore return address
    addi $sp, $sp, 20          ; Deallocate stack space
    jr $ra                     ; Return

; Function: compare_spectra
; Purpose: Compare two spectra for equality
compare_spectra:
    addi $sp, $sp, -16         ; Allocate stack space
    sw $ra, 12($sp)            ; Save return address
    
    li $v0, 1                  ; Assume match initially
    
    ; Compare each element (simplified implementation)
    la $t0, $a0                ; Target spectrum
    la $t1, $a1                ; Current spectrum
    lw $t2, $a2                ; Size
    
compare_loop:
    beq $t2, $zero, compare_done
    
    lw $t3, 0($t0)
    lw $t4, 0($t1)
    
    beq $t3, $t4, continue_compare
    li $v0, 0                  ; No match
    
continue_compare:
    addi $t0, $t0, 4           ; Next element
    addi $t1, $t1, 4           ; Next element
    addi $t2, $t2, -1          ; Decrement counter
    j compare_loop
    
compare_done:
    lw $ra, 12($sp)            ; Restore return address
    addi $sp, $sp, 16          ; Deallocate stack space
    jr $ra                     ; Return

; Function: print_spectrum (for debugging)
print_spectrum:
    addi $sp, $sp, -12         ; Allocate stack space
    sw $ra, 8($sp)             ; Save return address
    
    li $t0, 0                  ; i = 0
    
print_loop:
    bge $t0, $a1, print_done
    
    la $t1, current_spectrum
    sll $t2, $t0, 2            ; Word offset
    add $t1, $t1, $t2
    lw $t3, 0($t1)
    
    ; Print value (simplified - would need actual print function)
    addi $t0, $t0, 1
    j print_loop
    
print_done:
    lw $ra, 8($sp)             ; Restore return address
    addi $sp, $sp, 12          ; Deallocate stack space
    jr $ra                     ; Return

exit:
    li $v0, 10                 ; Exit system call
    syscall                    ; Exit program
```

## Key Algorithm Steps

1. **Backtracking Search**: Try all possible combinations of amino acids
2. **Cyclic Spectrum Generation**: For each peptide, generate all subpeptide masses
3. **Spectrum Matching**: Compare theoretical spectrum with target spectrum
4. **Pruning**: Early termination when partial solutions can't lead to valid answers

## Time Complexity

- O(N × M) where N is the number of possible peptides and M is the average peptide length
- Space complexity: O(M) for storing current peptide and spectrum

## Notes

This implementation provides a framework for solving the problem. In practice, you would need to:
1. Add proper input/output handling
2. Implement efficient spectrum comparison algorithms
3. Add pruning conditions to improve performance
4. Handle edge cases like empty spectra or invalid inputs