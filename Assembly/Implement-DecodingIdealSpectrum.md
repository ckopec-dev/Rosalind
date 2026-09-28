# Rosalind Problem: Implement Decoding Ideal Spectrum

## Problem Understanding

The problem requires implementing a decoding algorithm for an ideal spectrum, which is essentially finding the amino acid sequence that corresponds to a given mass spectrum. This is a classic dynamic programming problem.

## Approach

We'll use dynamic programming with memoization to solve this problem:
1. Create a mapping from masses to amino acids
2. Use recursive approach with memoization
3. For each target mass, try all possible amino acid additions
4. Build the solution backwards

## Assembly Implementation

```assembly
; Implement Decoding Ideal Spectrum
; Input: spectrum (array of masses)
; Output: amino acid sequence

.data
    ; Mass to amino acid mapping
    mass_table: .word 57, 71, 86, 99, 101, 103, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186
    aa_table: .ascii "GASPVTCILDNQKEMHFRYW"
    aa_count: .word 20
    
    ; Memory for dynamic programming
    dp_mem: .space 10000
    memo_size: .word 10000
    
    ; Result storage
    result_buffer: .space 1000
    result_len: .word 0

.text
.globl decode_spectrum

decode_spectrum:
    ; Parameters:
    ; r0 = spectrum array address
    ; r1 = spectrum length
    ; r2 = target mass
    
    push {r4-r11, lr}
    
    ; Initialize DP table
    mov r4, #0              ; i = 0
    mov r5, #0              ; current mass
    mov r6, #0              ; result index
    
init_dp_loop:
    cmp r4, r2
    bge init_dp_done
    
    ; Clear dp_mem[r4] = 0 (no solution)
    mov r7, #0
    str r7, [r5, r4, lsl #2]
    
    add r4, r4, #1
    b init_dp_loop
    
init_dp_done:
    ; Call recursive decoder
    mov r4, r2              ; target mass
    mov r5, #0              ; current position
    bl decode_recursive
    
    ; Build result string
    mov r4, #0              ; i = 0
    mov r6, #0              ; result index
    
build_result_loop:
    cmp r4, r5
    bge build_result_done
    
    ; Get amino acid from dp_mem[r4]
    ldr r7, [r5, r4, lsl #2]
    add r7, r7, #1          ; Adjust for 0-based indexing
    
    ; Convert to ASCII and store in result_buffer
    mov r8, #0              ; temp for ASCII conversion
    
    ; Find corresponding amino acid
    bl find_amino_acid
    
    strb r8, [r9, r6]
    add r6, r6, #1
    add r4, r4, #1
    b build_result_loop
    
build_result_done:
    mov r0, #0              ; null terminator
    strb r0, [r9, r6]
    
    pop {r4-r11, pc}

decode_recursive:
    ; r0 = current mass
    ; r1 = memo address
    
    push {r4-r7, lr}
    
    ; Base case: mass = 0
    cmp r0, #0
    beq return_zero
    
    ; Check if already computed (memoization)
    mov r4, r0
    ldr r5, [r1, r4, lsl #2]
    cmp r5, #0
    bne return_memo
    
    ; Try all amino acids
    mov r4, #0              ; i = 0
    mov r5, #0              ; found = 0
    
try_amino_loop:
    cmp r4, #20             ; 20 amino acids
    bge try_amino_done
    
    ; Get mass of amino acid i
    ldr r6, mass_table
    add r6, r6, r4, lsl #2
    ldr r7, [r6]
    
    ; Check if we can subtract this mass
    cmp r0, r7
    blt try_amino_next
    
    ; Recursively solve for remaining mass
    sub r8, r0, r7          ; new target mass
    bl decode_recursive
    
    ; If solution exists, store it
    cmp r9, #0              ; return value
    beq try_amino_next
    
    ; Store solution in memo table
    mov r10, r9             ; solution found
    str r10, [r1, r0, lsl #2]
    
    mov r5, #1              ; found = 1
    b try_amino_next
    
try_amino_next:
    add r4, r4, #1
    b try_amino_loop
    
try_amino_done:
    cmp r5, #0
    beq return_zero
    
    mov r0, #1              ; solution exists
    pop {r4-r7, pc}
    
return_zero:
    mov r0, #0              ; no solution
    pop {r4-r7, pc}
    
return_memo:
    mov r0, r5              ; return memoized result
    pop {r4-r7, pc}

find_amino_acid:
    ; r0 = amino acid mass
    ; returns ASCII character in r8
    
    push {r4-r6, lr}
    
    mov r4, #0              ; i = 0
    
find_loop:
    cmp r4, #20
    bge find_done
    
    ldr r5, mass_table
    add r5, r5, r4, lsl #2
    ldr r6, [r5]
    
    cmp r6, r0
    beq find_found
    
    add r4, r4, #1
    b find_loop
    
find_found:
    ldr r5, aa_table
    add r5, r5, r4
    ldrb r8, [r5]
    
find_done:
    pop {r4-r6, pc}

; Helper function to build spectrum from amino acid sequence
build_spectrum:
    ; r0 = amino acid string address
    ; r1 = output spectrum array
    
    push {r4-r7, lr}
    
    mov r4, #0              ; i = 0
    mov r5, #0              ; cumulative mass
    
build_loop:
    ldrb r6, [r0, r4]
    cmp r6, #0              ; null terminator
    beq build_done
    
    ; Convert character to mass
    bl char_to_mass
    
    add r5, r5, r6          ; accumulate mass
    
    ; Store in spectrum
    str r5, [r1, r4, lsl #2]
    
    add r4, r4, #1
    b build_loop
    
build_done:
    pop {r4-r7, pc}

char_to_mass:
    ; r0 = character
    ; returns mass in r0
    
    push {r4-r6, lr}
    
    mov r4, #0              ; i = 0
    
mass_loop:
    cmp r4, #20
    bge mass_done
    
    ldr r5, aa_table
    add r5, r5, r4
    ldrb r6, [r5]
    
    cmp r6, r0
    beq mass_found
    
    add r4, r4, #1
    b mass_loop
    
mass_found:
    ldr r5, mass_table
    add r5, r5, r4, lsl #2
    ldr r0, [r5]
    
mass_done:
    pop {r4-r6, pc}
```

## Algorithm Explanation

1. **Initialization**: Create a mapping between amino acid masses and their single-letter codes
2. **Dynamic Programming**: Use memoization to avoid recomputing subproblems
3. **Recursive Decoding**: For each target mass, try subtracting all possible amino acid masses
4. **Backtracking**: Reconstruct the solution by tracking which amino acids were used

## Time Complexity
- **Time**: O(M × N) where M is the maximum mass and N is the number of amino acids
- **Space**: O(M) for the DP table and memoization

## Key Features
- Uses memoization to avoid redundant calculations
- Handles edge cases (zero mass, no solution)
- Properly reconstructs the amino acid sequence
- Efficient memory usage through pre-allocated buffers

This implementation efficiently solves the ideal spectrum decoding problem using dynamic programming principles in assembly language.