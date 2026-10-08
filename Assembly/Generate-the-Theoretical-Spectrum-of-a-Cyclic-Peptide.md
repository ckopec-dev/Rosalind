# Rosalind Problem: Generate the Theoretical Spectrum of a Cyclic Peptide

## Problem Description
Given a cyclic peptide, generate its theoretical spectrum (the masses of all subpeptides).

## Solution in Assembly

```assembly
; Generate_theoretical_spectrum_of_a_cyclic_peptide
; Input: peptide sequence (string of amino acid symbols)
; Output: theoretical spectrum (sorted list of masses)

.data
    ; Amino acid masses (standard masses)
    aa_masses: .word 57, 71, 87, 97, 99, 101, 103, 113, 114, 115
               .word 128, 129, 131, 137, 147, 156, 163, 186
    
    ; Amino acid symbols
    aa_symbols: .ascii "GALMFWSPVTICYHNDKQER"
    
    ; Buffer for storing spectrum values
    spectrum_buffer: .space 1000
    
    ; Input peptide string (null terminated)
    peptide_input: .space 100
    
    ; Variables
    peptide_length: .long 0
    spectrum_count: .long 0

.text
.globl _start

_start:
    ; Initialize variables
    movl $0, %eax          ; spectrum_count = 0
    movl %eax, spectrum_count
    
    ; Get peptide length (assuming input is in peptide_input)
    call get_peptide_length
    movl %eax, peptide_length
    
    ; Generate all subpeptides for cyclic peptide
    call generate_cyclic_spectrum
    
    ; Sort spectrum values
    call sort_spectrum
    
    ; Print results
    call print_spectrum
    
    ; Exit program
    movl $1, %eax          ; sys_exit
    movl $0, %ebx          ; exit status
    int $0x80

; Function to get length of peptide string
get_peptide_length:
    pushl %ebp
    movl %esp, %ebp
    
    movl $peptide_input, %esi
    movl $0, %ecx          ; counter
    
length_loop:
    lodsb                  ; load byte from esi to al
    testb %al, %al         ; check if null terminator
    jz length_done
    incl %ecx              ; increment counter
    jmp length_loop
    
length_done:
    movl %ecx, %eax        ; return length
    popl %ebp
    ret

; Function to generate cyclic spectrum
generate_cyclic_spectrum:
    pushl %ebp
    movl %esp, %ebp
    
    movl peptide_length, %ecx
    movl $0, %edi          ; subpeptide_start_index
    
outer_loop:
    cmpb $0, %cl           ; if length = 0, exit
    jz outer_done
    
    ; Generate all subpeptides starting at index %edi
    call generate_subpeptides_at_index
    
    incl %edi              ; next start index
    decl %ecx              ; decrement remaining length
    jmp outer_loop
    
outer_done:
    popl %ebp
    ret

; Function to generate subpeptides starting at given index
generate_subpeptides_at_index:
    pushl %ebp
    movl %esp, %ebp
    
    movl peptide_length, %ecx
    movl $0, %esi          ; subpeptide_length
    
inner_loop:
    cmpb $0, %cl           ; if length = 0, exit
    jz inner_done
    
    ; Calculate mass of subpeptide from index %edi with length %esi+1
    call calculate_subpeptide_mass
    
    ; Store mass in spectrum_buffer
    movl %eax, spectrum_buffer(,%edx,4)  ; store at position edx*4
    
    incl %esi              ; increment subpeptide length
    decl %ecx              ; decrement remaining length
    jmp inner_loop
    
inner_done:
    popl %ebp
    ret

; Function to calculate mass of a subpeptide
calculate_subpeptide_mass:
    pushl %ebp
    movl %esp, %ebp
    
    ; Parameters: start_index in %edi, length in %esi
    movl $0, %eax          ; total mass = 0
    
    movl %esi, %ecx        ; copy length to ecx
    
subpeptide_loop:
    cmpb $0, %cl           ; if length = 0, exit
    jz subpeptide_done
    
    ; Get amino acid at current position (handle cyclic nature)
    call get_amino_acid_at_position
    movl %eax, %ebx        ; amino acid value
    
    ; Look up mass in aa_masses table
    movl aa_masses(,%ebx,4), %edx  ; get mass for this amino acid
    
    addl %edx, %eax        ; accumulate mass
    
    decb %cl               ; decrement length
    jmp subpeptide_loop
    
subpeptide_done:
    popl %ebp
    ret

; Function to get amino acid index at position (handle cyclic)
get_amino_acid_at_position:
    pushl %ebp
    movl %esp, %ebp
    
    ; Parameter: position in %edi
    movl $0, %eax          ; result = 0
    
    ; Convert position to array index with cyclic wraparound
    movl peptide_length, %ecx
    xorl %edx, %edx        ; clear edx for division
    divl %ecx              ; position / length (remainder in %edx)
    
    ; Now %edx contains the actual array index
    movb aa_symbols(,%edx,1), %al  ; get character
    
    ; Convert character to numeric index
    movl $0, %ebx          ; index = 0
    
char_loop:
    cmpb %al, aa_symbols(,%ebx,1)  ; compare with current symbol
    je char_found
    incl %ebx              ; next position
    jmp char_loop
    
char_found:
    movl %ebx, %eax        ; return index
    popl %ebp
    ret

; Function to sort spectrum (simple bubble sort)
sort_spectrum:
    pushl %ebp
    movl %esp, %ebp
    
    movl spectrum_count, %ecx
    
    ; Simple bubble sort
    movl $0, %esi          ; i = 0
    
bubble_outer:
    cmpb $0, %cl           ; if length <= 1, done
    jz bubble_done
    
    movl $0, %edi          ; j = 0
    
bubble_inner:
    ; Compare spectrum_buffer[j] with spectrum_buffer[j+1]
    movl spectrum_buffer(,%edi,4), %eax
    movl spectrum_buffer(4,%edi,4), %ebx
    
    cmpb %al, %bl          ; compare bytes
    jle bubble_continue    ; if already sorted, continue
    
    ; Swap elements
    movl %eax, %edx
    movl %ebx, spectrum_buffer(,%edi,4)
    movl %edx, spectrum_buffer(4,%edi,4)
    
bubble_continue:
    incl %edi              ; j++
    decl %ecx              ; decrement counter
    jmp bubble_inner
    
    decb %cl               ; decrement outer counter
    jmp bubble_outer
    
bubble_done:
    popl %ebp
    ret

; Function to print spectrum
print_spectrum:
    pushl %ebp
    movl %esp, %ebp
    
    movl spectrum_count, %ecx
    
print_loop:
    cmpb $0, %cl           ; if count = 0, exit
    jz print_done
    
    ; Print current value
    movl spectrum_buffer(,%edx,4), %eax
    call print_number
    
    ; Add space
    movl $32, %eax         ; space character
    call print_char
    
    decb %cl               ; decrement counter
    incl %edx              ; next element
    jmp print_loop
    
print_done:
    popl %ebp
    ret

; Helper function to print a number
print_number:
    pushl %ebp
    movl %esp, %ebp
    
    ; Simple implementation for demonstration
    ; In actual assembly, would need proper number-to-string conversion
    
    popl %ebp
    ret

; Helper function to print character
print_char:
    pushl %ebp
    movl %esp, %ebp
    
    ; Print character in %al
    ; This would use system calls in real implementation
    
    popl %ebp
    ret
```

## Key Concepts Implemented

1. **Cyclic Peptide Handling**: The algorithm accounts for the cyclic nature of the peptide by wrapping around the sequence when calculating subpeptides.

2. **Subpeptide Generation**: All possible contiguous subsequences are generated, including those that wrap around the end of the peptide.

3. **Mass Calculation**: Each amino acid is converted to its corresponding mass using a lookup table.

4. **Spectrum Sorting**: The resulting spectrum values are sorted in ascending order.

5. **Memory Management**: Uses buffers to store intermediate results and final spectrum.

## Time Complexity
- O(n³) where n is the length of the peptide
- Generating all subpeptides: O(n²) 
- Calculating masses: O(n) per subpeptide

## Space Complexity
- O(n²) for storing the spectrum values

This assembly implementation provides a framework for solving the cyclic peptide theoretical spectrum problem, though in practice it would need to be adapted to the specific assembly dialect and system calls available.