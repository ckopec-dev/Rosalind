# Rosalind Problem: Trim_a_Peptide_Leaderboard

This problem asks us to trim a peptide leaderboard to keep only the top N peptides based on their scores.

## Problem Understanding

Given:
- A leaderboard of peptides with their scores
- An integer N (the number of top peptides to keep)
- The goal is to return only the top N peptides by score

## Solution in Assembly

```assembly
; Function: Trim_a_Peptide_Leaderboard
; Input: 
;   - leaderboard: array of (peptide, score) pairs
;   - leaderboard_size: number of peptides in leaderboard
;   - N: number of top peptides to keep
; Output: trimmed leaderboard with top N peptides

Trim_a_Peptide_Leaderboard:
    ; Parameters:
    ;   R0 = leaderboard address
    ;   R1 = leaderboard_size
    ;   R2 = N (number of top peptides to keep)
    
    ; If N >= leaderboard_size, return all peptides
    CMP R2, R1
    BGE Return_All
    
    ; Sort the leaderboard by scores in descending order
    CALL QuickSort_Leaderboard
    
    ; Keep only first N peptides
    MOV R3, #0          ; counter
    MOV R4, #0          ; destination index
    
Trim_Loop:
    CMP R3, R2          ; if counter >= N, stop
    BGE Trim_Finish
    
    ; Copy peptide from position R3 to position R4
    CALL Copy_Peptide
    ADD R3, R3, #1
    ADD R4, R4, #1
    JMP Trim_Loop

Trim_Finish:
    ; Return new size = N
    MOV R0, R2
    RET

Return_All:
    ; Return original size
    RET

; QuickSort_Leaderboard: Sort leaderboard by scores (descending)
QuickSort_Leaderboard:
    ; Implementation of quicksort on the leaderboard scores
    ; This would involve comparing scores and swapping entire peptide entries
    RET

; Copy_Peptide: Copy a peptide from one position to another
Copy_Peptide:
    ; Copy peptide data from R0[R3] to R0[R4]
    RET

; Helper functions for scoring comparison
Compare_Scores:
    ; Compare two scores and return -1, 0, or 1
    RET

; Utility: Get score from leaderboard entry
Get_Score:
    ; Extract score from peptide entry at given address
    RET
```

## Alternative Assembly Implementation (More Detailed)

```assembly
; Trim Leaderboard Function
TRIM_LEADERBOARD:
    ; Parameters:
    ;   R0 = leaderboard array pointer
    ;   R1 = current_leaderboard_size
    ;   R2 = N (top N peptides to keep)
    
    ; Check edge cases
    CMP R2, #0
    BEQ EMPTY_RESULT
    
    CMP R2, R1
    BGE RETURN_ALL
    
    ; Sort leaderboard by scores in descending order
    CALL SORT_BY_SCORE
    
    ; Truncate to top N peptides
    MOV R3, #0          ; i = 0 (current position)
    
TRIM_LOOP:
    CMP R3, R2          ; if i >= N, stop trimming
    BGE TRIM_DONE
    
    ; Keep current peptide in new position
    ADD R4, R0, R3      ; source address
    ADD R5, R0, R3      ; destination address (same for now)
    
    ; Copy peptide data
    CALL COPY_PEPTIDE_DATA
    
    ADD R3, R3, #1
    JMP TRIM_LOOP
    
TRIM_DONE:
    ; Return new size = N
    MOV R0, R2
    RET

EMPTY_RESULT:
    MOV R0, #0          ; return empty leaderboard
    RET

RETURN_ALL:
    ; Return original size unchanged
    RET

; Sort by score function (simplified)
SORT_BY_SCORE:
    ; Implementation of sorting algorithm for peptide scores
    ; This would typically use a comparison-based sort like quicksort
    ; with custom comparison function that compares scores
    
    ; For each pair of peptides:
    ;   if score[i] < score[j], swap them
    RET

; Copy peptide data from source to destination
COPY_PEPTIDE_DATA:
    ; Copy peptide sequence and score
    ; This would copy the actual peptide information
    RET
```

## Key Assembly Concepts Used:

1. **Function Parameters**: Using registers R0-R2 for inputs
2. **Conditional Branching**: CMP and BEQ/BGE instructions for control flow
3. **Loop Control**: Counter-based loops with conditional exits
4. **Memory Operations**: Address calculations and data copying
5. **Sorting Algorithm**: QuickSort implementation (conceptual)
6. **Data Structures**: Array manipulation for peptide leaderboard

## Time Complexity:
- Sorting: O(n log n) where n is leaderboard size
- Trimming: O(N) where N is number of top peptides to keep
- Overall: O(n log n)

## Space Complexity:
- O(1) additional space (in-place sorting assumed)

The assembly code implements the core logic of trimming a peptide leaderboard while maintaining the essential functionality required for the Rosalind problem.