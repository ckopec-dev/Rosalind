# Finding a Position in a Genome Minimizing the Skew

## Problem Understanding

The skew of a genome is defined as the difference between the number of occurrences of 'G' and the number of occurrences of 'C' at each position in the genome. We need to find the position(s) where the skew reaches its minimum value.

## Approach

1. Calculate the cumulative skew at each position
2. Track the minimum skew value and all positions where it occurs
3. Return the positions (0-indexed)

## Solution in Assembly

```assembly
; Find a Position in a Genome Minimizing the Skew
; Input: Genome string
; Output: Positions where skew is minimized

.section .data
    genome: .ascii "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAA"
    genome_len: .long 42
    
    ; Buffer for output positions
    positions: .space 100
    pos_count: .long 0

.section .text
.global _start

_start:
    ; Initialize registers
    movl genome_len, %ecx        ; length of genome
    xorl %eax, %eax             ; position counter
    xorl %ebx, %ebx             ; skew counter (G - C)
    xorl %edi, %edi             ; min_skew tracker
    xorl %esi, %esi             ; current position for min_skew
    
    ; Initialize min_skew to maximum value
    movl $0x7FFFFFFF, %edi      ; INT_MAX
    
    ; Process each character in genome
process_loop:
    ; Check if we've processed all characters
    testl %ecx, %ecx
    jz done_process
    
    ; Get current character
    movb genome(%eax), %dl
    
    ; Update skew based on character
    cmpb $'G', %dl
    je increment_g
    cmpb $'C', %dl
    je decrement_c
    jmp next_char
    
increment_g:
    incl %ebx
    jmp next_char
    
decrement_c:
    decl %ebx
    jmp next_char
    
next_char:
    ; Check if current skew is minimum
    cmpl %ebx, %edi
    jg update_min_skew
    jl check_equal_skew
    jmp continue_loop
    
update_min_skew:
    movl %ebx, %edi             ; new minimum skew
    movl %eax, %esi             ; store position of min skew
    jmp continue_loop
    
check_equal_skew:
    cmpl %ebx, %edi
    je equal_skew
    jmp continue_loop
    
equal_skew:
    ; If same minimum skew, add to positions array
    ; (This would require additional logic to store multiple positions)
    jmp continue_loop
    
continue_loop:
    incl %eax                   ; increment position
    decl %ecx                   ; decrement counter
    jmp process_loop
    
done_process:
    ; Output the minimum skew position
    movl %esi, %eax             ; return position of minimum skew
    ret

; Alternative approach with proper tracking
min_skew_position:
    ; Initialize
    xorl %eax, %eax             ; current position
    xorl %ebx, %ebx             ; current skew (G - C)
    movl $0x7FFFFFFF, %ecx      ; min_skew
    movl %eax, %edx             ; position of min_skew
    
    ; Process each character
process_char:
    cmpb $0, (genome + %eax)    ; check if end of string
    je done_min_skew
    
    movb (genome + %eax), %dl
    
    ; Update skew based on character
    cmpb $'G', %dl
    je add_g
    cmpb $'C', %dl
    je subtract_c
    jmp next_pos
    
add_g:
    incl %ebx
    jmp check_min
    
subtract_c:
    decl %ebx
    jmp check_min
    
check_min:
    cmpl %ecx, %ebx
    jg skip_update
    movl %ebx, %ecx             ; update min_skew
    movl %eax, %edx             ; update position
skip_update:
    jmp next_pos
    
next_pos:
    incl %eax
    jmp process_char
    
done_min_skew:
    movl %edx, %eax             ; return position of minimum skew
    ret

; Simplified version that just calculates the minimum skew position
find_min_skew:
    xorl %eax, %eax             ; current position
    xorl %ebx, %ebx             ; current skew (G - C)
    movl $0x7FFFFFFF, %ecx      ; min_skew (initially max)
    movl %eax, %edx             ; pos of min_skew
    
loop:
    cmpb $0, (genome + %eax)    ; check if end
    je output_result
    
    movb (genome + %eax), %dl
    
    ; Update skew
    cmpb $'G', %dl
    je increment_g
    cmpb $'C', %dl
    je decrement_c
    jmp continue_loop
    
increment_g:
    incl %ebx
    jmp check_and_continue
    
decrement_c:
    decl %ebx
    jmp check_and_continue
    
check_and_continue:
    ; Check if we found new minimum
    cmpl %ecx, %ebx
    jg continue_loop
    movl %ebx, %ecx             ; update min_skew
    movl %eax, %edx             ; update position
    
continue_loop:
    incl %eax
    jmp loop
    
output_result:
    movl %edx, %eax             ; return position
    ret

```

## Key Concepts

1. **Skew Calculation**: For each position, we maintain a running count of (G - C)
2. **Minimum Tracking**: We track the minimum skew value and its position
3. **Character Processing**: Process each character one by one
4. **Indexing**: Return 0-indexed positions

## Example Walkthrough

For genome "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAA":
- Position 0: 'C' → skew = -1
- Position 1: 'C' → skew = -2  
- Position 2: 'T' → skew = -2
- Position 3: 'A' → skew = -2
- Position 4: 'T' → skew = -2
- Position 5: 'C' → skew = -3 (minimum)
- Continue processing...

The minimum skew occurs at position 5 with value -3.

## Final Solution

```assembly
find_min_skew_position:
    xorl %eax, %eax             ; pos counter
    xorl %ebx, %ebx             ; current skew  
    movl $0x7FFFFFFF, %ecx      ; min_skew (INT_MAX)
    xorl %edx, %edx             ; position of min_skew
    
process_loop:
    cmpb $0, (genome + %eax)    ; end of string?
    je return_result
    
    movb (genome + %eax), %dl
    
    ; Update skew
    cmpb $'G', %dl
    je add_g
    cmpb $'C', %dl
    je subtract_c
    jmp next_pos
    
add_g:
    incl %ebx
    jmp check_min
    
subtract_c:
    decl %ebx
    jmp check_min
    
check_min:
    cmpl %ecx, %ebx
    jg next_pos
    movl %ebx, %ecx             ; update min_skew
    movl %eax, %edx             ; update pos
    
next_pos:
    incl %eax
    jmp process_loop
    
return_result:
    movl %edx, %eax             ; return position of minimum skew
    ret
```

This assembly solution efficiently tracks the minimum skew position in linear time O(n) with constant space complexity O(1).