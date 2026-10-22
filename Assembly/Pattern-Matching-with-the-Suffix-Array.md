# Rosalind Problem: Pattern Matching with the Suffix Array (Assembly Solution)

## Problem Understanding

Given a string and a pattern, we need to find all occurrences of the pattern in the string using suffix array information.

## Approach

We'll implement a binary search approach on the suffix array to find pattern matches efficiently.

## Assembly Implementation

```assembly
.section .data
    # Input string and pattern (example)
    text: .ascii "ACGTACGTACGT"
    text_len: .long 12
    
    pattern: .ascii "CGT"
    pattern_len: .long 3
    
    newline: .ascii "\n"
    
    # Suffix array (precomputed for example)
    suffix_array: .long 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
    sa_len: .long 12

.section .text
    .global _start

_start:
    # Initialize registers
    movl text_len, %ecx          # Length of text
    movl pattern_len, %edx       # Length of pattern
    
    # Call pattern matching function
    call pattern_match_suffix_array
    
    # Exit program
    movl $1, %eax                # sys_exit
    movl $0, %ebx                # exit status
    int $0x80

# Function: pattern_match_suffix_array
# Input: text, pattern, suffix_array (global)
# Output: positions where pattern occurs in text
pattern_match_suffix_array:
    pushl %ebp
    movl %esp, %ebp
    
    # Initialize search boundaries
    movl $0, %esi                # low = 0
    movl sa_len, %edi            # high = suffix array length
    decl %edi                    # high = sa_len - 1
    
    # Binary search for lower bound
    call binary_search_lower_bound
    
    # Store lower bound result
    movl %eax, %ebx              # lower_bound position
    
    # Binary search for upper bound
    call binary_search_upper_bound
    
    # Store upper bound result  
    movl %eax, %ecx              # upper_bound position
    
    # Print results between bounds
    call print_matches
    
    popl %ebp
    ret

# Binary search for lower bound of pattern in suffix array
binary_search_lower_bound:
    pushl %ebp
    movl %esp, %ebp
    
    movl %esi, %eax              # low
    movl %edi, %ebx              # high
    
binary_search_loop_lower:
    cmpl %ebx, %eax              # if low > high
    jg binary_search_end_lower
    
    # Calculate mid = (low + high) / 2
    addl %ebx, %eax
    shrl $1, %eax                # right shift by 1 = divide by 2
    
    # Get suffix at position mid
    movl suffix_array(,%eax,4), %ecx  # get suffix start position
    movl text, %edx              # text address
    
    # Compare pattern with suffix starting at position %ecx
    pushl %eax                   # save mid
    call compare_suffix_with_pattern
    popl %eax                    # restore mid
    
    cmpl $0, %eax                # if comparison result == 0
    je binary_search_end_lower   # found match, continue search to left
    
    jl binary_search_lower_bound_less  # pattern < suffix, search left
    jmp binary_search_lower_bound_greater  # pattern > suffix, search right

binary_search_lower_bound_less:
    # pattern < suffix, search left
    decl %eax                    # mid - 1
    movl %eax, %ebx              # high = mid - 1
    jmp binary_search_loop_lower

binary_search_lower_bound_greater:
    # pattern > suffix, search right
    incl %eax                    # mid + 1
    movl %eax, %esi              # low = mid + 1
    jmp binary_search_loop_lower

binary_search_end_lower:
    # Return the final position
    movl %eax, %eax              # return value
    
    popl %ebp
    ret

# Binary search for upper bound of pattern in suffix array
binary_search_upper_bound:
    pushl %ebp
    movl %esp, %ebp
    
    movl %esi, %eax              # low
    movl %edi, %ebx              # high
    
binary_search_loop_upper:
    cmpl %ebx, %eax              # if low > high
    jg binary_search_end_upper
    
    # Calculate mid = (low + high) / 2
    addl %ebx, %eax
    shrl $1, %eax                # right shift by 1 = divide by 2
    
    # Get suffix at position mid
    movl suffix_array(,%eax,4), %ecx  # get suffix start position
    movl text, %edx              # text address
    
    # Compare pattern with suffix starting at position %ecx
    pushl %eax                   # save mid
    call compare_suffix_with_pattern
    popl %eax                    # restore mid
    
    cmpl $0, %eax                # if comparison result == 0
    je binary_search_loop_upper  # found match, continue search to right
    
    jl binary_search_upper_bound_less  # pattern < suffix, search left
    jmp binary_search_upper_bound_greater  # pattern > suffix, search right

binary_search_upper_bound_less:
    # pattern < suffix, search left
    decl %eax                    # mid - 1
    movl %eax, %ebx              # high = mid - 1
    jmp binary_search_loop_upper

binary_search_upper_bound_greater:
    # pattern > suffix, search right
    incl %eax                    # mid + 1
    movl %eax, %esi              # low = mid + 1
    jmp binary_search_loop_upper

binary_search_end_upper:
    # Return the final position
    movl %eax, %eax              # return value
    
    popl %ebp
    ret

# Compare suffix with pattern
# Input: text address in %edx, suffix start position in %ecx
# Output: comparison result in %eax (0=equal, <0, >0)
compare_suffix_with_pattern:
    pushl %ebp
    movl %esp, %ebp
    
    # Copy pattern to stack for comparison
    movl pattern, %edi           # pattern address
    movl %ecx, %esi              # suffix start position in text
    
    movl pattern_len, %ecx       # pattern length
    
compare_loop:
    cmpl $0, %ecx                # if length = 0
    je compare_equal
    
    # Get character from suffix and pattern
    movb text(,%esi,1), %al      # get suffix char
    movb pattern(,%edi,1), %bl   # get pattern char
    
    cmpl %bl, %al                # compare characters
    jne compare_not_equal        # if not equal, return difference
    
    incl %esi                    # increment suffix position
    incl %edi                    # increment pattern position
    decl %ecx                    # decrement counter
    
    jmp compare_loop
    
compare_equal:
    movl $0, %eax                # strings are equal
    jmp compare_end
    
compare_not_equal:
    movl %al, %eax               # return difference
    subl %bl, %eax               # calculate difference
    
compare_end:
    popl %ebp
    ret

# Print matches between lower and upper bounds
print_matches:
    pushl %ebp
    movl %esp, %ebp
    
    # Simple loop to print positions
    movl %ebx, %esi              # start position (lower bound)
    movl %ecx, %edi              # end position (upper bound)
    
print_loop:
    cmpl %edi, %esi              # if start >= end
    jge print_end
    
    # Print position (simplified - in real implementation would format properly)
    pushl %esi                   # position to print
    call print_number            # function to print number
    popl %esi
    
    # Print newline
    pushl $newline
    call print_char
    popl %esi
    
    incl %esi                    # next position
    jmp print_loop
    
print_end:
    popl %ebp
    ret

# Simple print number function (simplified)
print_number:
    pushl %ebp
    movl %esp, %ebp
    
    # Convert number to string and print
    # Implementation would involve division and digit conversion
    
    popl %ebp
    ret

# Simple print character function (simplified)
print_char:
    pushl %ebp
    movl %esp, %ebp
    
    # Write character to stdout
    # This is a simplified version - real implementation would use system calls
    
    popl %ebp
    ret
```

## Key Assembly Concepts Used

1. **Binary Search Implementation**: Efficient pattern matching using binary search on precomputed suffix array
2. **String Comparison**: Character-by-character comparison of text and pattern
3. **Memory Management**: Proper use of stack and registers for function calls
4. **Control Flow**: Conditional jumps and loops for search logic
5. **System Calls**: Basic exit functionality using Linux system calls

## Time Complexity
- O(log n + k) where n is the text length and k is the number of matches
- Binary search on suffix array takes O(log n)
- Printing matches takes O(k)

This assembly implementation demonstrates how to perform pattern matching using suffix arrays in a low-level language while maintaining the algorithmic efficiency of the approach.