# Rosalind Problem: Compute the Hamming Distance Between Two Strings (Assembly Solution)

## Problem Understanding

The Hamming distance between two strings of equal length is the number of positions at which the corresponding symbols are different.

## Assembly Solution

```assembly
.data
    str1: .ascii "GAGCCTACTAACGGGAT\n"
    str2: .ascii "CATCGGTATCGGCCCTA\n"
    len: .long 17
    result: .long 0

.text
    .global _start

_start:
    # Initialize pointers and counter
    movl str1, %esi        # pointer to first string
    movl str2, %edi        # pointer to second string
    movl len, %ecx         # length counter
    xorl %eax, %eax        # result counter (hamming distance)
    
hamming_loop:
    # Check if we've reached the end
    testl %ecx, %ecx
    jz hamming_done
    
    # Compare characters at current position
    movb (%esi), %al       # load character from first string
    movb (%edi), %bl       # load character from second string
    
    # If characters are different, increment counter
    cmpb %bl, %al
    jne increment_counter
    jmp next_position
    
increment_counter:
    incl %eax              # increment hamming distance
    
next_position:
    incl %esi              # move to next character in first string
    inc %edi               # move to next character in second string
    decl %ecx              # decrement counter
    jmp hamming_loop
    
hamming_done:
    # Store result
    movl %eax, result
    
    # Exit program
    movl $1, %eax          # sys_exit
    movl $0, %ebx          # exit status
    int $0x80
```

## Alternative Implementation (More Efficient)

```assembly
.data
    str1: .ascii "GAGCCTACTAACGGGAT\n"
    str2: .ascii "CATCGGTATCGGCCCTA\n"
    len: .long 17

.text
    .global _start

_start:
    movl str1, %esi        # pointer to first string
    movl str2, %edi        # pointer to second string
    movl len, %ecx         # length counter
    xorl %eax, %eax        # hamming distance counter
    
hamming_loop:
    testl %ecx, %ecx
    jz hamming_done
    
    # Compare bytes directly
    cmpb (%esi), (%edi)
    jne increment_counter
    
increment_counter:
    incl %eax              # increment if different
    
next_char:
    incl %esi
    incl %edi
    decl %ecx
    jmp hamming_loop
    
hamming_done:
    # Result is now in %eax
    # For demonstration, we'll just exit with the result
    movl $1, %eax          # sys_exit
    movl %eax, %ebx        # use result as exit code
    int $0x80
```

## Explanation

This assembly program computes the Hamming distance between two DNA strings:

1. **Data Section**: 
   - Two DNA strings of equal length (17 characters)
   - Length counter
   - Result storage location

2. **Algorithm**:
   - Initialize pointers to both strings and length counter
   - Loop through each character position
   - Compare characters at current position using `cmpb`
   - If different, increment the hamming distance counter
   - Continue until all positions are checked

3. **Key Instructions**:
   - `movb`: Move bytes between memory and registers
   - `cmpb`: Compare bytes (sets flags)
   - `jne`: Jump if not equal (increment counter)
   - `incl`: Increment long integer
   - `testl`: Test if counter is zero

4. **Result**: The Hamming distance (number of different positions) is stored in the result variable and returned as exit code.

For the given example strings, the Hamming distance would be 7.