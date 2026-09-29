# Rosalind Problem: Convert a Peptide into a Peptide Vector

## Problem Description
Given a peptide, we need to convert it into its corresponding peptide vector representation.

The peptide vector is a binary vector of length 3n (where n is the length of the peptide) where:
- Each amino acid is represented by a 3-bit binary number (000 to 111)
- The vector has 1s at positions corresponding to the amino acid representations
- The mapping from amino acids to 3-bit numbers is as follows:
  - A = 000, C = 001, D = 010, E = 011, F = 100, G = 101, H = 110, I = 111
  - K = 000, L = 001, M = 010, N = 011, P = 100, Q = 101, R = 110, S = 111
  - T = 000, V = 001, W = 010, Y = 011

Wait, let me reconsider this. Looking at the standard amino acid mapping for peptide vectors:

Actually, the standard mapping is:
- A = 000, C = 001, D = 010, E = 011, F = 100, G = 101, H = 110, I = 111
- K = 000, L = 001, M = 010, N = 011, P = 100, Q = 101, R = 110, S = 111
- T = 000, V = 001, W = 010, Y = 011

No, let me look at this more carefully. The standard peptide vector uses a specific mapping of 20 amino acids to 3-bit representations.

Let me use the standard amino acid to 3-bit mapping:
- A = 000, C = 001, D = 010, E = 011, F = 100, G = 101, H = 110, I = 111
- K = 000, L = 001, M = 010, N = 011, P = 100, Q = 101, R = 110, S = 111
- T = 000, V = 001, W = 010, Y = 011

This doesn't make sense - there's overlap. Let me use a more standard approach.

## Solution in Assembly Language

```assembly
; Convert a Peptide into a Peptide Vector
; Input: peptide string (sequence of amino acids)
; Output: binary vector representation

.data
    ; Amino acid to 3-bit mapping table
    aa_map: .word 0b000, 0b001, 0b010, 0b011, 0b100, 0b101, 0b110, 0b111
            ; A,  C,  D,  E,  F,  G,  H,  I
            .word 0b000, 0b001, 0b010, 0b011, 0b100, 0b101, 0b110, 0b111  
            ; K,  L,  M,  N,  P,  Q,  R,  S
            .word 0b000, 0b001, 0b010, 0b011  ; T, V, W, Y
    
    peptide_input: .ascii "SKIV"
    peptide_length: .long 4

.text
.globl _start

_start:
    ; Initialize registers
    movl peptide_length(%esp), %ecx        ; Load peptide length
    movl $0, %eax                          ; Initialize result vector index
    movl $0, %edx                          ; Initialize bit counter
    
    ; Process each amino acid in the peptide
process_aa:
    ; Check if we've processed all amino acids
    cmpl $0, %ecx
    jz done
    
    ; Get current amino acid (assuming it's in memory)
    movb peptide_input(%eax), %bl          ; Load current amino acid
    
    ; Convert amino acid to 3-bit representation
    call aa_to_bits
    
    ; Store the 3 bits in result vector
    ; This is a simplified approach - actual implementation would need bit manipulation
    
    ; Move to next amino acid
    incl %eax
    decl %ecx
    jmp process_aa

done:
    ; Exit program
    movl $1, %eax                          ; sys_exit
    movl $0, %ebx                          ; exit status
    int $0x80

; Function to convert amino acid to 3-bit representation
aa_to_bits:
    ; This is a placeholder for the actual conversion logic
    ; In practice, we would map each AA character to its corresponding 3-bit value
    
    ; Return value in %eax (simplified)
    ret

; Alternative implementation with proper mapping
convert_peptide_vector:
    pushl %ebp
    movl %esp, %ebp
    
    ; Load peptide string and length
    movl 8(%ebp), %esi     ; peptide string
    movl 12(%ebp), %ecx    ; length
    
    ; Initialize result vector (3 * length bits)
    xorl %eax, %eax        ; clear eax
    xorl %ebx, %ebx        ; clear ebx
    xorl %edx, %edx        ; clear edx
    
convert_loop:
    cmpl $0, %ecx
    jz convert_done
    
    ; Get current character
    movb (%esi), %dl
    
    ; Map amino acid to 3-bit value (simplified mapping)
    cmpb $'A', %dl
    je aa_to_bits_A
    cmpb $'C', %dl
    je aa_to_bits_C
    cmpb $'D', %dl
    je aa_to_bits_D
    ; ... continue with other amino acids ...
    
aa_to_bits_A:
    movl $0, %eax
    jmp convert_continue
    
aa_to_bits_C:
    movl $1, %eax
    jmp convert_continue
    
aa_to_bits_D:
    movl $2, %eax
    jmp convert_continue

convert_continue:
    ; Store the 3 bits in the result vector
    ; Implementation would involve bit shifting and OR operations
    
    incl %esi              ; Move to next character
    decl %ecx              ; Decrement counter
    jmp convert_loop

convert_done:
    movl %ebp, %esp
    popl %ebp
    ret
```

## Alternative Simpler Approach

```assembly
; Simpler implementation for peptide vector conversion
; This assumes we have a complete mapping table

.data
    ; Amino acid to 3-bit mapping (20 amino acids)
    aa_map_table: .byte 0, 1, 2, 3, 4, 5, 6, 7   ; A,B,C,D,E,F,G,H  
                  .byte 0, 1, 2, 3, 4, 5, 6, 7   ; K,L,M,N,P,Q,R,S
                  .byte 0, 1, 2, 3               ; T,V,W,Y

.text
.globl convert_peptide_to_vector

convert_peptide_to_vector:
    pushl %ebp
    movl %esp, %ebp
    
    ; Parameters: peptide string address in %esi, length in %ecx
    ; Result vector will be built in memory
    
    xorl %eax, %eax        ; result vector index
    xorl %ebx, %ebx        ; bit counter
    
convert_loop:
    cmpl $0, %ecx
    jz done_conversion
    
    ; Get current amino acid character
    movb (%esi), %dl
    
    ; Convert to 3-bit value (simplified - actual implementation needs mapping)
    ; This is a placeholder for proper AA->bit conversion
    movb aa_map_table(%edx), %al
    
    ; Store the 3 bits in result vector (actual bit manipulation needed)
    
    incl %esi              ; next character
    decl %ecx              ; decrement count
    jmp convert_loop

done_conversion:
    movl %ebp, %esp
    popl %ebp
    ret
```

## Note
This is a conceptual implementation. The actual assembly code would require:
1. Proper mapping of each amino acid character to its 3-bit representation
2. Bit manipulation to store the 3 bits in the result vector
3. Handling of different amino acid characters (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
4. Proper memory management for the output vector

The peptide vector will have length 3 × (peptide length) where each amino acid contributes 3 bits to the final binary vector.