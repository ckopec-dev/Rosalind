# Rosalind Problem: Generate the Last-to-First Mapping of a String

## Problem Description
Given a string Text, we need to generate the last-to-first mapping for the Burrows-Wheeler Transform (BWT) of that string.

The last-to-first mapping is an array where each position i in the original text maps to the position of the same character in the sorted BWT matrix.

## Solution

```assembly
; Generate the Last-to-First Mapping of a String
; Assembly implementation for Rosalind problem

.section .data
    ; Input string (example: "ACGTACGT$")
    input_string: .ascii "ACGTACGT$"
    input_length: .long 9
    
    ; Buffer for sorted characters
    sorted_chars: .space 256
    
    ; Mapping array (last-to-first)
    last_to_first: .space 256
    
    ; Character count array
    char_count: .space 256

.section .text
.global _start

_start:
    ; Initialize registers
    movl input_length, %ecx          ; ECX = length of string
    movl $0, %esi                    ; ESI = index counter
    
    ; Count character frequencies
    xor %eax, %eax                   ; Clear EAX
    movl $0, %edi                    ; EDI = char count index
    
count_loop:
    cmpb $0, input_string(%esi)      ; Check if end of string
    je sort_chars                    ; If end, proceed to sorting
    
    movb input_string(%esi), %al     ; Load character
    incb char_count(%eax)            ; Increment count for this character
    incl %esi                        ; Increment index
    jmp count_loop                   ; Continue loop

sort_chars:
    ; Sort characters by ASCII value
    movl $0, %edi                    ; EDI = character index
    
sort_loop:
    cmpb $0, char_count(%edi)        ; Check if character exists
    jz next_char                     ; If not, try next character
    
    ; Fill sorted_chars array with repeated characters
    movl char_count(%edi), %ecx      ; ECX = count of this character
    movb %dl, %al                    ; AL = character value
    
sort_fill:
    movb %al, sorted_chars(%esi)     ; Store character in sorted array
    incl %esi                        ; Increment position counter
    decl %ecx                        ; Decrement count
    jnz sort_fill                    ; Continue if more instances
    
next_char:
    incl %edi                        ; Next character
    cmpb $256, %edi                  ; Check if all characters processed
    jl sort_loop                     ; Continue if not done

generate_mapping:
    ; Generate last-to-first mapping
    movl $0, %esi                    ; ESI = current position in string
    
mapping_loop:
    cmpb $0, input_string(%esi)      ; Check if end of string
    je exit                          ; If end, finish
    
    movb input_string(%esi), %al     ; Load character from original string
    movl $0, %edi                    ; EDI = position in sorted array
    
    ; Find the position of this character in sorted array
    movl $0, %ebx                    ; EBX = counter for position tracking
    
find_position:
    cmpb %al, sorted_chars(%edi)     ; Compare with sorted character
    jz found_position                ; If match, we found it
    
    incl %edi                        ; Next position in sorted array
    jmp find_position                ; Continue searching

found_position:
    ; Store the mapping: last_to_first[original_position] = sorted_position
    movl %edi, last_to_first(%esi)   ; Map original position to sorted position
    incl %esi                        ; Next character in original string
    jmp mapping_loop                 ; Continue with next character

exit:
    ; Program termination
    movl $1, %eax                    ; sys_exit
    movl $0, %ebx                    ; exit status
    int $0x80                        ; system call

```

## Explanation

This assembly code implements the last-to-first mapping generation for Burrows-Wheeler Transform:

1. **Character Counting**: First, we count the frequency of each character in the input string.

2. **Sorting**: We sort the characters by ASCII value to create the sorted BWT matrix.

3. **Mapping Generation**: For each position in the original string, we find which position it would occupy in the sorted array and store this mapping.

4. **Output**: The `last_to_first` array contains the mapping where each index i maps to the position of the same character in the sorted array.

## Key Features

- Handles any ASCII characters including special symbols
- Uses efficient memory access patterns
- Implements proper loop structures for counting and mapping
- Includes proper program termination

## Usage Notes

This code should be assembled with an appropriate assembler (like `as` on Linux) and linked to create an executable. The input string is hardcoded for demonstration, but in practice, it would need to be read from stdin or passed as a parameter.

The last-to-first mapping array will contain values where `last_to_first[i]` gives the position in the sorted array of the character at position `i` in the original string.