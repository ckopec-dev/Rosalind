# Rosalind Problem: Implement PatternToNumber in Assembly

## Problem Understanding

PatternToNumber converts a DNA pattern (string of nucleotides) into its lexicographically ordered index in the list of all possible k-mers sorted lexicographically.

For example, with k=3:
- Pattern "AAA" → Index 0
- Pattern "AAT" → Index 1
- Pattern "ATA" → Index 4
- Pattern "TTT" → Index 63

## Assembly Solution

```assembly
.data
    pattern BYTE "ACGT", 0        ; Input pattern (null-terminated)
    k DWORD 4                     ; Length of pattern
    
.code
PatternToNumber PROC
    ; Input: pointer to null-terminated DNA string in ECX
    ; Output: DWORD result in EAX
    ; Uses: EAX, EBX, ECX, EDX, ESI
    
    push ebp
    mov ebp, esp
    push ebx
    push ecx
    push edx
    push esi
    
    mov esi, ecx                ; esi = pointer to pattern string
    xor eax, eax                ; result = 0
    xor ebx, ebx                ; temp variable
    
    ; Get length of pattern
    mov ecx, k                  ; load k (length)
    
    ; Process each character from left to right
    mov edx, 0                  ; power counter
    
    ; Loop through each position in the pattern
    pattern_loop:
        cmp ecx, 0              ; if we've processed all characters
        je pattern_done
        
        ; Get current character
        mov bl, [esi]           ; load current character
        cmp bl, 0               ; check for null terminator
        je pattern_done
        
        ; Convert character to numeric value (A=0, C=1, G=2, T=3)
        xor edx, edx            ; clear edx
        mov dl, bl              ; move character to dl
        cmp dl, 'A'             ; check if it's A
        je is_A
        cmp dl, 'C'             ; check if it's C
        je is_C
        cmp dl, 'G'             ; check if it's G
        je is_G
        cmp dl, 'T'             ; check if it's T
        je is_T
        
        ; Invalid character - should not happen in valid input
        jmp invalid_char
        
    is_A:
        xor edx, edx            ; edx = 0
        jmp convert_done
        
    is_C:
        mov edx, 1              ; edx = 1
        jmp convert_done
        
    is_G:
        mov edx, 2              ; edx = 2
        jmp convert_done
        
    is_T:
        mov edx, 3              ; edx = 3
        
    convert_done:
        ; Multiply current result by 4 and add new digit
        mov ebx, eax            ; save old result
        shl eax, 2              ; multiply by 4 (eax = result * 4)
        add eax, edx            ; add the new digit
        
        ; Move to next character
        inc esi                 ; advance pointer
        dec ecx                 ; decrease counter
        jmp pattern_loop
        
    pattern_done:
        ; Result is in EAX
        pop esi
        pop edx
        pop ecx
        pop ebx
        mov esp, ebp
        pop ebp
        ret
        
    invalid_char:
        ; Handle error case - return 0 or appropriate error value
        xor eax, eax
        jmp pattern_done

PatternToNumber ENDP

; Alternative implementation using lookup table
.data
    char_to_num BYTE 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ; Initialize lookup table for ASCII values
    ; This would normally be set up in a more complex way

.code
PatternToNumberOptimized PROC
    ; Optimized version using lookup table approach
    push ebp
    mov ebp, esp
    push ebx
    push ecx
    push edx
    
    mov esi, ecx                ; esi = pointer to pattern string
    xor eax, eax                ; result = 0
    
    ; Process each character from left to right
    mov ecx, k                  ; load k (length)
    
    process_loop:
        cmp ecx, 0              ; if we've processed all characters
        je done_processing
        
        ; Get current character
        mov bl, [esi]           ; load current character
        
        ; Convert character to numeric value using lookup table
        ; For simplicity, we'll use the direct approach:
        ; A=0, C=1, G=2, T=3
        cmp bl, 'A'
        je add_A
        cmp bl, 'C'
        je add_C
        cmp bl, 'G'
        je add_G
        cmp bl, 'T'
        je add_T
        
        ; Invalid character - error handling
        jmp invalid_input
        
    add_A:
        mov edx, 0
        jmp convert_done2
        
    add_C:
        mov edx, 1
        jmp convert_done2
        
    add_G:
        mov edx, 2
        jmp convert_done2
        
    add_T:
        mov edx, 3
        
    convert_done2:
        ; Multiply result by 4 and add new digit
        shl eax, 2              ; multiply result by 4 (eax = result * 4)
        add eax, edx            ; add current digit
        
        inc esi                 ; move to next character
        dec ecx                 ; decrease counter
        jmp process_loop
        
    done_processing:
        pop edx
        pop ecx
        pop ebx
        mov esp, ebp
        pop ebp
        ret
        
    invalid_input:
        xor eax, eax            ; return 0 for error case
        jmp done_processing

PatternToNumberOptimized ENDP

END
```

## Key Algorithm Concepts

1. **Base-4 Number System**: DNA uses 4 nucleotides (A,C,G,T) so this is essentially converting from base-4 to decimal.

2. **Lexicographic Ordering**: The algorithm processes characters from left to right, treating each position as a digit in base-4.

3. **Positional Value**: Each character contributes to the final result based on its position and value:
   - First character: contributes `value × 4^(k-1)`
   - Second character: contributes `value × 4^(k-2)`
   - And so on...

## Example Walkthrough

For pattern "ACGT" (k=4):
1. A = 0, position 3: 0 × 4³ = 0
2. C = 1, position 2: 1 × 4² = 16  
3. G = 2, position 1: 2 × 4¹ = 8
4. T = 3, position 0: 3 × 4⁰ = 3

Total = 0 + 16 + 8 + 3 = 27

## Assembly Considerations

- Uses registers efficiently (EAX for result, ECX for loop counter)
- Handles null-terminated strings
- Implements character-to-number conversion
- Performs bit shifting for multiplication by powers of 4
- Includes error handling for invalid input characters