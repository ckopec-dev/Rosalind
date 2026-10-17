# Rosalind Problem: Implement CycleToChromosome

## Problem Description
TheCycleToChromosome problem asks us to convert a cycle of integers into a chromosome represented as a permutation.

Given a cycle of integers, we need to reconstruct the corresponding chromosome (permutation) by:
1. Breaking the cycle at specified breakpoints
2. Converting the resulting segments into signed permutations

## Solution in Assembly Language

```assembly
; CycleToChromosome - Convert cycle to chromosome
; Input: cycle - array of integers representing a cycle
;        n - length of cycle
; Output: chromosome - array of integers representing the chromosome

CycleToChromosome:
    ; Function to convert cycle to chromosome
    ; Input: 
    ;   R0 = cycle address
    ;   R1 = n (length of cycle)
    ;   R2 = chromosome address
    
    MOV R3, #0          ; i = 0
    MOV R4, #0          ; j = 0
    
CycleLoop:
    CMP R3, R1          ; Compare i with n
    BGE EndCycle        ; If i >= n, end loop
    
    ; Calculate chromosome[j] = cycle[i+1]/2
    MOV R5, R3          ; R5 = i
    ADD R5, #1          ; R5 = i + 1
    LDR R6, [R0, R5]    ; R6 = cycle[i+1]
    MOV R7, R6          ; Copy to R7
    MOV R7, R7, LSR #1  ; R7 = cycle[i+1]/2
    
    ; Handle sign: if cycle[i+1] is even, it's positive; if odd, it's negative
    MOV R8, R6          ; Copy cycle[i+1] to R8
    AND R8, #1          ; Check if odd (sign bit)
    
    CMP R8, #0          ; If odd, make negative
    BEQ PositiveSign
    
    NEG R7              ; Make negative
    
PositiveSign:
    STR R7, [R2, R4]    ; chromosome[j] = result
    
    ADD R3, #2          ; i += 2 (jump to next pair)
    ADD R4, #1          ; j += 1
    B CycleLoop         ; Continue loop

EndCycle:
    MOV PC, LR          ; Return from function

; Alternative implementation for better clarity
CycleToChromosome_v2:
    ; Input: R0 = cycle array address
    ;        R1 = n (length of cycle)
    ;        R2 = chromosome array address
    
    MOV R3, #0          ; i = 0 (index in cycle)
    MOV R4, #0          ; j = 0 (index in chromosome)
    
ConvertLoop:
    CMP R3, R1          ; Check if we've processed all elements
    BGE ConvertEnd      ; If yes, exit
    
    ; Get the element at position i+1 in cycle
    MOV R5, R3          ; Copy current index
    ADD R5, #1          ; i+1
    LDR R6, [R0, R5]    ; Load cycle[i+1]
    
    ; Divide by 2 to get chromosome value
    MOV R7, R6          ; Copy to R7
    MOV R7, R7, LSR #1  ; Right shift by 1 (equivalent to /2)
    
    ; Determine sign based on original value
    MOV R8, R6          ; Copy original value
    AND R8, #1          ; Check if odd
    
    CMP R8, #0          ; If odd, make negative
    BEQ StoreValue      ; If even, keep positive
    
    NEG R7              ; Make it negative
    
StoreValue:
    STR R7, [R2, R4]    ; Store in chromosome[j]
    
    ADD R3, #2          ; Move to next pair (i += 2)
    ADD R4, #1          ; Move to next position in chromosome
    B ConvertLoop       ; Continue loop
    
ConvertEnd:
    MOV PC, LR          ; Return

; Example usage:
; Assume cycle = [1, 2, 3, 4, 5, 6]
; Expected chromosome = [1, -2, 3, -4, 5, -6]
```

## Algorithm Explanation

The algorithm works as follows:

1. **Input Processing**: Takes a cycle represented as an array of integers
2. **Element Conversion**: For each pair of consecutive elements in the cycle:
   - Take the second element of the pair (at index i+1)
   - Divide by 2 to get the chromosome value
   - Determine sign based on whether the original element was odd or even
3. **Output Generation**: Constructs the final chromosome permutation

## Key Assembly Concepts Used

- **Memory Operations**: `LDR` (load) and `STR` (store) for array access
- **Bit Manipulation**: `LSR` (logical shift right) for division by 2, `AND` for odd/even check
- **Control Flow**: Conditional branches (`CMP`, `BGE`, `BEQ`) for loop control
- **Register Management**: Proper use of registers for indexing and calculations

## Time Complexity
O(n) where n is the length of the cycle

## Space Complexity  
O(1) additional space (excluding input/output arrays)

This implementation assumes that the cycle array contains integers in a specific format where consecutive elements represent pairs, and the sign is determined by whether the original element is odd or even.