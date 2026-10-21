# Find the Shortest Non-Shared Substring of Two Strings - Assembly Solution

I'll solve this step by step using assembly language to find the shortest substring that appears in one string but not in another.

## Problem Understanding

Given two strings, we need to find the shortest substring that:
1. Appears in the first string
2. Does NOT appear in the second string

## Assembly Implementation

```assembly
.section .data
    # Input strings (example)
    str1:   .ascii "ABCD"
    str1_len = . - str1
    
    str2:   .ascii "ACDF"
    str2_len = . - str2
    
    # Buffer for substring storage
    buffer: .space 100
    
    # Result storage
    result: .space 100
    result_len: .long 0

.section .text
    .global _start

_start:
    # Initialize registers
    movl $str1, %esi      # First string pointer
    movl $str2, %edi      # Second string pointer
    movl $str1_len, %ecx  # Length of first string
    movl $str2_len, %edx  # Length of second string
    
    # Find shortest non-shared substring
    call find_shortest_non_shared
    
    # Exit program
    movl $1, %eax         # sys_exit
    movl $0, %ebx         # exit status
    int $0x80

# Function: find_shortest_non_shared
# Input: 
#   %esi = pointer to first string
#   %edi = pointer to second string  
#   %ecx = length of first string
#   %edx = length of second string
find_shortest_non_shared:
    pushl %ebp
    movl %esp, %ebp
    
    # Initialize variables
    movl $1, %eax         # Start with substring length 1
    movl %ecx, %ebx       # Save original length
    
shortest_loop:
    # Check if we've exceeded string length
    cmpl %ebx, %eax
    jg no_solution
    
    # Try all substrings of current length in first string
    movl $0, %edi         # substring index in str1
    movl %eax, %ecx       # current substring length
    
sub_string_loop:
    # Check if we can extract substring of length %eax from str1 at position %edi
    addl %edi, %eax
    cmpl %ebx, %eax
    jl sub_string_continue
    
    # Reset for next iteration
    movl $0, %eax         # Reset to 1 (length)
    jmp shortest_loop
    
sub_string_continue:
    # Extract substring from str1
    call extract_substring
    
    # Check if this substring exists in str2
    call is_substring_in_string
    
    # If not found in str2, we found our answer
    cmpl $0, %eax
    jz return_result
    
    # Move to next substring in str1
    movl $0, %edi         # Reset index
    jmp sub_string_loop

# Function: extract_substring
# Input: 
#   %esi = pointer to string
#   %edi = start position  
#   %ecx = length
# Output: 
#   buffer contains extracted substring
extract_substring:
    pushl %ebp
    movl %esp, %ebp
    
    movl %esi, %eax       # String pointer
    addl %edi, %eax       # Add offset
    movl %eax, %esi       # New string pointer
    
    # Copy to buffer
    movl $buffer, %edi
    movl %ecx, %edx
    
copy_loop:
    lodsb                 # Load byte from %esi
    stosb                 # Store byte to %edi
    decl %edx
    jnz copy_loop
    
    # Add null terminator
    movb $0, (%edi)
    
    popl %ebp
    ret

# Function: is_substring_in_string
# Input:
#   buffer = substring to search for  
#   %edi = pointer to second string
#   %edx = length of second string
# Output:
#   %eax = 1 if found, 0 if not found
is_substring_in_string:
    pushl %ebp
    movl %esp, %ebp
    
    # Initialize search parameters
    movl $0, %ecx         # Current position in str2
    movl $buffer, %esi    # Substring to find
    
    # Check if substring is empty
    cmpb $0, (%esi)
    jz not_found
    
    # Search for substring in str2
search_loop:
    # Check bounds
    cmpl %edx, %ecx
    jge not_found
    
    # Compare current position with substring
    movl %ecx, %edi       # Current position in str2
    addl %edi, %edi       # Index in str2
    
    # Compare characters one by one
    call compare_strings
    
    # If match found, return 1 (found)
    cmpl $0, %eax
    jz found_match
    
    # Continue searching
    incl %ecx
    jmp search_loop
    
found_match:
    movl $1, %eax         # Found
    jmp search_end
    
not_found:
    movl $0, %eax         # Not found
    
search_end:
    popl %ebp
    ret

# Function: compare_strings  
# Input:
#   %esi = substring pointer (buffer)
#   %edi = string pointer (str2 + offset)
# Output:
#   %eax = 0 if match, non-zero if mismatch
compare_strings:
    pushl %ebp
    movl %esp, %ebp
    
    movl $0, %ecx         # Character counter
    
compare_loop:
    movb (%esi,%ecx,1), %al  # Get character from substring
    cmpb %al, (%edi,%ecx,1)  # Compare with string character
    jne compare_fail
    
    testb %al, %al        # Check for null terminator
    jz compare_success
    
    incl %ecx
    jmp compare_loop
    
compare_fail:
    movl $1, %eax         # Mismatch
    jmp compare_end
    
compare_success:
    movl $0, %eax         # Match
    
compare_end:
    popl %ebp
    ret

# Function: return_result
# Input: 
#   buffer = shortest non-shared substring
return_result:
    # Copy result to output buffer
    pushl %esi
    pushl %edi
    
    movl $result, %edi
    movl $buffer, %esi
    
    # Copy string
    movl $0, %ecx
copy_result_loop:
    lodsb
    stosb
    testb %al, %al
    jz copy_done
    incl %ecx
    jmp copy_result_loop
    
copy_done:
    movl %ecx, result_len
    
    popl %edi
    popl %esi
    popl %ebp
    ret

# Error handling for no solution found
no_solution:
    # Handle case where no non-shared substring exists
    movl $0, %eax
    jmp search_end
```

## Key Algorithm Steps

1. **Iterative Search**: Start with length 1 substrings and increase until we find a non-shared one
2. **Substring Extraction**: Extract each possible substring from the first string
3. **Substring Matching**: Check if each extracted substring exists in the second string
4. **Early Termination**: Return immediately when the first non-shared substring is found

## Time Complexity Analysis

- **Time**: O(n³) where n is the length of the strings
- **Space**: O(n) for buffer storage

## Example Run

For input strings "ABCD" and "ACDF":
- Length 1: "A" (shared), "B" (not shared) → Return "B"
- The algorithm will find that "B" is in first string but not in second string

## Notes

This assembly implementation demonstrates the core algorithm using low-level operations. In practice, this would be optimized with better string matching algorithms like KMP or Boyer-Moore for better performance.