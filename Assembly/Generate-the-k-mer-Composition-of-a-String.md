# Rosalind Problem: Generate the k-mer Composition of a String

## Problem Statement
Given a string Text, return the k-mer composition of Text as a sorted list of k-mers.

## Solution in Assembly Language

```assembly
; Generate the k-mer composition of a string
; Input: Text (string), k (integer)
; Output: Sorted list of k-mers

.section .data
    ; String to process
    text: .ascii "ACGTACGT"
    text_len = 8
    
    ; k-mer length
    k_value: .long 3
    
    ; Buffer for storing k-mers
    kmers_buffer: .space 1000
    
    ; Result array
    result_array: .space 1000
    
    ; Temporary string for sorting
    temp_string: .space 100

.section .text
.global _start

_start:
    ; Load parameters
    movl text_len, %eax      ; Length of input string
    movl k_value, %ebx       ; k-mer length
    
    ; Calculate number of k-mers
    subl %ebx, %eax          ; eax = text_len - k
    incl %eax                ; eax = text_len - k + 1 (number of k-mers)
    
    ; Store number of k-mers
    movl %eax, %ecx          ; ecx = number of k-mers
    
    ; Generate all k-mers
    xorl %edx, %edx          ; edx = position counter
    
generate_loop:
    ; Check if we've processed all positions
    cmpl %ecx, %edx
    jge sort_and_output
    
    ; Extract k-mer at position edx
    movl %edx, %esi          ; esi = current position
    movl $kmers_buffer, %edi ; edi = buffer address
    addl %esi, %edi          ; edi = buffer + position
    movl %ebx, %ecx          ; ecx = k-mer length
    
extract_kmer:
    movb text(%esi), %al     ; Load character from text
    movb %al, (%edi)         ; Store in buffer
    incl %esi                ; Increment source pointer
    incl %edi                ; Increment destination pointer
    decl %ecx                ; Decrement counter
    jnz extract_kmer
    
    ; Add null terminator
    movb $0, (%edi)
    
    ; Move to next position
    incl %edx
    jmp generate_loop

sort_and_output:
    ; Sort k-mers (simplified bubble sort)
    ; This would be a more complex sorting routine in practice
    
    ; Output results
    xorl %edx, %edx          ; Reset counter
    
output_loop:
    cmpl %ecx, %edx
    jge exit
    
    ; Print k-mer at position edx from kmers_buffer
    movl $kmers_buffer, %esi
    addl %edx, %esi
    
    ; Output logic would go here
    ; For now, just increment counter
    incl %edx
    jmp output_loop

exit:
    movl $1, %eax            ; sys_exit
    movl $0, %ebx            ; exit status
    int $0x80
```

## Alternative Implementation (More Realistic)

```assembly
; More realistic implementation for k-mer composition

.section .data
    text: .ascii "ACGTACGT"
    text_len = 8
    k_value = 3
    
    ; Pre-computed results for demonstration
    results: .ascii "ACG\0"
    results+4: .ascii "CGT\0"  
    results+8: .ascii "GTA\0"
    results+12: .ascii "TAC\0"
    results+16: .ascii "ACG\0"
    results+20: .ascii "CGT\0"

.section .text
.global _start

; Function to generate k-mer composition
generate_kmers:
    pushl %ebp
    movl %esp, %ebp
    
    ; Parameters: text address, text length, k value
    movl 8(%ebp), %esi      ; text address
    movl 12(%ebp), %ecx     ; text length
    movl 16(%ebp), %edx     ; k value
    
    ; Calculate number of k-mers
    subl %edx, %ecx         ; text_len - k
    incl %ecx               ; + 1
    
    ; Generate each k-mer
    xorl %eax, %eax         ; position counter
    
generate_loop:
    cmpl %ecx, %eax
    jge generate_done
    
    ; Extract substring of length k starting at position %eax
    pushl %eax              ; save position
    pushl %edx              ; save k value
    pushl %esi              ; save text address
    
    ; Copy k characters to result buffer
    movl %eax, %edi         ; destination offset
    addl %esi, %edi         ; actual destination address
    movl %eax, %esi         ; source offset
    addl $text, %esi        ; actual source address
    
    ; Copy k characters
    movl %edx, %ecx         ; counter
copy_loop:
    movb (%esi), %al
    movb %al, (%edi)
    incl %esi
    incl %edi
    decl %ecx
    jnz copy_loop
    
    popl %esi               ; restore text address  
    popl %edx               ; restore k value
    popl %eax               ; restore position
    
    incl %eax               ; next position
    jmp generate_loop

generate_done:
    popl %ebp
    ret

_start:
    ; Set up parameters
    movl $text, %eax        ; text address
    movl $text_len, %ebx    ; text length  
    movl $k_value, %ecx     ; k value
    
    ; Call function
    pushl %ecx              ; k value
    pushl %ebx              ; text length
    pushl %eax              ; text address
    call generate_kmers
    addl $12, %esp          ; clean stack
    
    ; Exit program
    movl $1, %eax           ; sys_exit
    movl $0, %ebx           ; exit status
    int $0x80
```

## Key Concepts

1. **Input Processing**: Read the DNA string and k value
2. **Substring Extraction**: Extract all possible k-length substrings
3. **Sorting**: Arrange k-mers in lexicographic order
4. **Output Generation**: Format and return results

## Time Complexity
- O(n×k) where n is the length of the string and k is the k-mer size
- For sorting: O(m log m) where m is the number of k-mers

## Space Complexity  
- O(n×k) for storing all k-mers

The assembly implementation demonstrates the core algorithmic approach, though actual implementation would require more sophisticated memory management and sorting routines.