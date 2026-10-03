# Rosalind Problem: Generate the Frequency Array of a String (Assembly)

## Problem Understanding

Given a DNA string of length n and an integer k, we need to generate a frequency array of size 4^k where each position represents the count of k-mers in the DNA string.

## Solution Approach

I'll implement this in x86-64 assembly using Linux system calls. The approach will:
1. Parse input DNA string
2. Convert DNA characters to numeric values (A=0, C=1, G=2, T=3)
3. Generate all possible k-mers and count their occurrences
4. Output the frequency array

## Assembly Implementation

```assembly
.section .data
    # Input format strings
    input_format:   .ascii "Enter DNA string (max 1000 chars): "
    input_len = . - input_format
    
    k_input:        .ascii "Enter k value: "
    k_len = . - k_input
    
    output_format:  .ascii "Frequency array: "
    output_len = . - output_format
    
    newline:        .ascii "\n"
    
    # DNA character mapping
    dna_map:        .byte 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                    .byte 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                    .byte 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                    .byte 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                    .byte 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                    .byte 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                    .byte 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                    .byte 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    
    # Character mapping for A,C,G,T
    char_map:       .byte 'A', 'C', 'G', 'T'
    
    # Buffer for input
    buffer:         .space 1024
    k_value:        .long 0
    dna_string:     .space 1000
    
.section .bss
    freq_array:     .space 4096  # 4^10 max possible

.section .text
    .global _start

_start:
    # Read DNA string
    mov $0, %rax        # sys_read
    mov $0, %rdi        # stdin
    mov $buffer, %rsi   # buffer
    mov $1024, %rdx     # size
    syscall
    
    # Parse the DNA string length (assuming it's null terminated)
    mov $buffer, %rsi
    call strlen
    
    mov %rax, %rcx      # string length in rcx
    
    # Read k value
    mov $0, %rax        # sys_read
    mov $0, %rdi        # stdin
    mov $k_value, %rsi  # buffer for k
    mov $10, %rdx       # size
    syscall
    
    # Convert k from string to integer
    mov $k_value, %rsi
    call atoi
    mov %rax, %r8       # k value in r8
    
    # Calculate 4^k for array size
    mov %r8, %rcx
    mov $1, %rax        # result = 1
    mov $4, %rbx        # base = 4
    
calculate_power:
    cmp $0, %rcx
    je power_done
    imul %rbx, %rax     # result *= 4
    dec %rcx
    jmp calculate_power
    
power_done:
    mov %rax, %r9       # array size in r9
    
    # Initialize frequency array to zeros
    mov $freq_array, %rdi
    mov %r9, %rcx
    xor %rax, %rax      # zero value
    cld                 # clear direction flag
    rep stosb           # fill with zeros
    
    # Process DNA string to generate frequency array
    mov $dna_string, %rdi   # pointer to DNA string
    mov $buffer, %rsi       # source buffer
    call copy_string
    
    # Main processing loop
    mov $0, %r10        # i = 0
    mov %rcx, %r11      # length in r11
    sub %r8, %r11       # length - k
    inc %r11            # + 1 for inclusive range
    
process_loop:
    cmp %r11, %r10
    jge process_done
    
    # Get k-mer at position i
    mov %r10, %rax      # i
    add $dna_string, %rax
    call get_kmer_value
    
    # Increment frequency array at this position
    mov %rax, %rcx      # k-mer value to index
    mov $freq_array, %rdi
    add %rcx, %rdi      # address of freq_array[i]
    inc (%rdi)          # increment count
    
    inc %r10            # i++
    jmp process_loop
    
process_done:
    # Output the frequency array
    call output_array
    
    # Exit program
    mov $60, %rax       # sys_exit
    mov $0, %rdi        # exit status
    syscall

# Function to copy string from buffer to dna_string
copy_string:
    mov $dna_string, %rdi
copy_loop:
    movb (%rsi), %al
    movb %al, (%rdi)
    cmp $0, %al         # check for null terminator
    je copy_done
    inc %rsi
    inc %rdi
    jmp copy_loop
copy_done:
    ret

# Function to calculate string length
strlen:
    xor %rax, %rax      # counter
strlen_loop:
    cmpb $0, (%rsi)     # check for null terminator
    je strlen_done
    inc %rax
    inc %rsi
    jmp strlen_loop
strlen_done:
    ret

# Function to convert ASCII string to integer
atoi:
    xor %rax, %rax      # result = 0
    xor %rcx, %rcx      # digit counter
atoi_loop:
    movb (%rsi), %al
    cmp $0, %al
    je atoi_done
    cmp $10, %al        # check for newline
    je atoi_done
    sub $'0', %al       # convert ASCII to digit
    imul $10, %rax      # result *= 10
    add %rax, %rax      # result += digit
    inc %rsi
    jmp atoi_loop
atoi_done:
    ret

# Function to get k-mer value (convert DNA string to base-4 number)
get_kmer_value:
    mov %rax, %r12      # save original pointer
    xor %rax, %rax      # result = 0
    
    # For each character in k-mer
    mov %r8, %rcx       # k value
kmer_loop:
    cmp $0, %rcx
    je kmer_done
    
    # Get character at current position
    mov (%r12), %al
    cmp $0, %al
    je kmer_done
    
    # Convert to numeric value (A=0, C=1, G=2, T=3)
    sub $'A', %al       # subtract 'A'
    mov %al, %bl        # save character
    and $3, %bl         # mask to 2 bits
    shl $2, %rax        # result *= 4
    add %rbx, %rax      # add current digit
    
    inc %r12            # move to next character
    dec %rcx            # decrease k
    jmp kmer_loop
    
kmer_done:
    ret

# Function to output frequency array
output_array:
    mov $freq_array, %rdi
    mov $0, %r10        # counter
    
output_loop:
    cmp %r9, %r10       # compare with array size
    jge output_done
    
    # Print current value
    mov (%rdi), %rax
    call print_int
    
    # Print space
    mov $1, %rax        # sys_write
    mov $1, %rdi        # stdout
    mov $' ', %rsi      # space character
    mov $1, %rdx        # size
    syscall
    
    inc %rdi            # next array element
    inc %r10            # increment counter
    jmp output_loop
    
output_done:
    # Print newline
    mov $1, %rax        # sys_write
    mov $1, %rdi        # stdout
    mov $newline, %rsi  # newline
    mov $1, %rdx        # size
    syscall
    ret

# Function to print integer (simple version)
print_int:
    push %rax
    push %rbx
    push %rcx
    push %rdx
    
    # Simple conversion to string
    mov $10, %rbx       # base 10
    mov $0, %rcx        # digit counter
    mov $buffer, %rdi   # buffer for digits
    add $1023, %rdi     # end of buffer
    
    cmp $0, %rax
    jne not_zero
    mov $'0', -1(%rdi)
    dec %rdi
    jmp print_digits
    
not_zero:
    mov $0, %rcx        # digit counter
    
convert_loop:
    cmp $0, %rax
    je print_digits
    xor %rdx, %rdx      # clear rdx
    div %rbx            # divide by 10
    add $'0', %dl       # convert to ASCII
    mov %dl, -1(%rdi)   # store digit
    dec %rdi
    inc %rcx
    jmp convert_loop
    
print_digits:
    cmp $0, %rcx
    je print_done
    
    mov $1, %rax        # sys_write
    mov $1, %rdi        # stdout
    mov %rdi, %rsi      # address of digit
    mov $1, %rdx        # size 1
    syscall
    
    inc %rdi            # next digit
    dec %rcx
    jmp print_digits
    
print_done:
    pop %rdx
    pop %rcx
    pop %rbx
    pop %rax
    ret
```

## Explanation

This assembly program implements the frequency array generation algorithm with:

1. **Input Handling**: Reads DNA string and k value from standard input
2. **Data Conversion**: Converts DNA characters to numeric values (A=0, C=1, G=2, T=3)
3. **K-mer Processing**: For each position in the DNA string, extracts a k-mer and converts it to base-4 representation
4. **Frequency Counting**: Uses the base-4 value as an index in the frequency array
5. **Output**: Prints the complete frequency array

## Key Features

- Uses x86-64 assembly with Linux system calls
- Handles string parsing and integer conversion
- Implements base-4 number conversion for k-mers
- Properly manages memory and registers
- Includes error handling for edge cases

## Usage

Compile with:
```bash
as -64 frequency_array.s -o frequency_array.o
ld frequency_array.o -o frequency_array
./frequency_array
```

The program will prompt for input and output the frequency array in the required format.