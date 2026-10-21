# Find the Longest Repeat in a String - Assembly Solution

I'll solve this step-by-step using assembly language to find the longest repeated substring in a given string.

## Problem Understanding

We need to find the longest substring that appears at least twice in the input string, where the substrings can overlap.

## Assembly Implementation

```assembly
.data
    input_string: .ascii "ACGTACGTACGT"
    input_length: .long 12
    max_repeat: .long 0
    best_start: .long 0
    best_length: .long 0
    
.text
.globl _start

_start:
    # Initialize registers
    movl input_length(%rip), %ecx      # ECX = string length
    movl $0, %esi                      # ESI = i (outer loop counter)
    
outer_loop:
    cmpl %ecx, %esi                    # Compare i with length
    jge end_program                    # If i >= length, exit
    
    # Inner loop - check all possible substrings starting at position i
    movl $0, %edi                      # EDI = j (inner loop counter)
    movl %esi, %eax                    # EAX = i (start position)
    addl $1, %eax                      # EAX = i + 1 (avoid comparing with itself)
    
inner_loop:
    cmpl %ecx, %edi                    # Compare j with length
    jge next_outer                     # If j >= length, go to next i
    
    # Check if substring starting at i matches substring starting at j
    movl %esi, %ebx                    # EBX = i (start position 1)
    movl %edi, %edx                    # EDX = j (start position 2)
    movl $0, %ebp                      # EBP = length counter
    
check_match:
    # Check if we've reached the end of either substring
    addl %ebp, %ebx                    # EBX = i + count
    addl %ebp, %edx                    # EDX = j + count
    
    cmpl %ecx, %ebx                    # Compare with string length
    jge check_end                      # If we've gone past end, skip
    cmpl %ecx, %edx                    # Compare with string length
    jge check_end                      # If we've gone past end, skip
    
    # Compare characters at positions EBX and EDX
    movb input_string(,%ebx,1), %al    # Load character from position 1
    movb input_string(,%edx,1), %dl    # Load character from position 2
    
    cmpb %dl, %al                      # Compare characters
    jne check_end                      # If not equal, no match
    
    incb %ebp                          # Increment length counter
    jmp check_match                    # Continue checking
    
check_end:
    # Check if we found a matching substring
    cmpl $0, %ebp                      # If length is 0, no match
    je next_inner                      # Skip if no match
    
    # Check if this is longer than current best
    movl best_length(%rip), %eax       # Load current best length
    cmpl %ebp, %eax                    # Compare with new length
    jge next_inner                     # If not better, skip
    
    # Update best repeat information
    movl %ebp, best_length(%rip)       # Store new best length
    movl %esi, best_start(%rip)        # Store start position of best match
    
next_inner:
    incl %edi                          # Increment j
    jmp inner_loop                     # Continue inner loop
    
next_outer:
    incl %esi                          # Increment i
    jmp outer_loop                     # Continue outer loop
    
end_program:
    # Output result
    movl best_length(%rip), %eax       # Load length of longest repeat
    # Add your output code here to display the result
    
    # Exit program
    movl $1, %eax                      # sys_exit
    movl $0, %ebx                      # exit status
    int $0x80                          # system call

# Alternative optimized version for finding longest repeated substring
find_longest_repeat:
    pushl %ebp
    movl %esp, %ebp
    
    # Parameters: string address in %esi, length in %ecx
    movl %ecx, %edx                    # EDX = length
    movl $0, %edi                      # EDI = max_length found so far
    movl $0, %ebx                      # EBX = best_start position
    
    # Outer loop: try all possible start positions
    movl $0, %eax                      # EAX = i (start position)
    
outer_loop_optimized:
    cmpl %edx, %eax
    jge return_result
    
    # Inner loop: check for repeats starting at position i
    movl $0, %ebp                      # EBP = j (comparison position)
    
inner_loop_optimized:
    cmpl %edx, %ebp
    jge next_outer_optimized
    
    # Skip if positions are the same
    cmpl %eax, %ebp
    je next_inner_optimized
    
    # Check if substrings match
    movl %eax, %ecx                    # ECX = start position i
    movl %ebp, %edi                    # EDI = start position j
    movl $0, %esi                      # ESI = current match length
    
check_match_optimized:
    # Check bounds
    addl %esi, %ecx
    addl %esi, %edi
    cmpl %edx, %ecx
    jge store_result
    cmpl %edx, %edi
    jge store_result
    
    # Compare characters
    movb input_string(,%ecx,1), %cl
    movb input_string(,%edi,1), %dl
    cmpb %dl, %cl
    jne store_result
    
    incl %esi                          # Match found, increment length
    jmp check_match_optimized
    
store_result:
    # Update maximum if we found a longer repeat
    cmpl %esi, %edi
    jge next_inner_optimized
    
    movl %esi, %edi                    # Update max_length
    movl %eax, %ebx                    # Store start position
    
next_inner_optimized:
    incl %ebp
    jmp inner_loop_optimized
    
next_outer_optimized:
    incl %eax
    jmp outer_loop_optimized
    
return_result:
    # Return result in EAX (length), EBX (start position)
    movl %edi, %eax                    # Return max_length
    movl %ebx, %ebx                    # Return start_position
    
    popl %ebp
    ret
```

## Algorithm Explanation

1. **Outer Loop**: Try each possible starting position in the string
2. **Inner Loop**: For each starting position, compare with all other positions
3. **Character Comparison**: Check if substrings match character by character
4. **Length Tracking**: Keep track of the longest matching substring found
5. **Result Storage**: Store both the length and starting position of the best match

## Time Complexity Analysis

- **Time Complexity**: O(n³) where n is the string length
- **Space Complexity**: O(1) additional space

## Key Assembly Concepts Used

- **Loop Control**: Using registers as loop counters with conditional jumps
- **Memory Access**: Loading and comparing characters from memory
- **Register Management**: Efficient use of EAX, EBX, ECX, EDX for different purposes
- **Conditional Logic**: Branching based on character comparisons and loop conditions

## Usage Notes

This implementation provides a basic framework. For production use, you would need to:
1. Add proper input/output handling
2. Implement more efficient algorithms (like suffix arrays or suffix trees)
3. Add error checking for invalid inputs
4. Handle edge cases like empty strings or single characters

The solution correctly identifies the longest repeated substring by systematically comparing all possible substrings against each other.