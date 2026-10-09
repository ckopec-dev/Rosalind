# Rosalind Problem: Compute the Number of Peptides of Given Total Mass

## Problem Understanding

We need to compute how many different peptides can have a given total mass, where each amino acid has a specific integer mass.

## Solution Approach

I'll use dynamic programming to count the number of ways to achieve each possible mass using the given amino acid masses.

## Assembly Implementation

```assembly
.section .data
# Amino acid masses (integer values)
.amino_masses:
    .long 57   # A
    .long 71   # C
    .long 87   # D
    .long 97   # E
    .long 99   # F
    .long 101  # G
    .long 103  # H
    .long 113  # I
    .long 114  # K
    .long 115  # L
    .long 128  # M
    .long 129  # N
    .long 131  # P
    .long 137  # Q
    .long 147  # R
    .long 156  # S
    .long 163  # T
    .long 186  # W
    .long 189  # Y
    .long 197  # V

# Number of amino acids
.amino_count: .long 20

.section .text
.globl _start

_start:
    # Initialize registers
    movl $20, %ecx          # Number of amino acids
    movl $1000, %edx        # Maximum possible mass (assuming input is <= 1000)
    
    # Allocate memory for dp array
    # dp[i] = number of ways to achieve mass i
    pushl %edx              # Save max_mass
    movl %edx, %eax
    addl $4, %eax           # +4 for array indexing
    pushl %eax              # Push size of dp array
    
    # Initialize dp array to zero
    xorl %esi, %esi         # i = 0
    movl %edx, %edi         # max_mass
    xorl %eax, %eax         # zero value for memset
    
memset_loop:
    movl %eax, (%esi)       # dp[i] = 0
    addl $4, %esi           # i++
    decl %edi
    jnz memset_loop
    
    # Base case: one way to achieve mass 0 (empty peptide)
    movl $1, (%esi)         # dp[0] = 1
    
    # Dynamic Programming
    movl $1, %esi           # i = 1 (start from mass 1)
    
dp_outer_loop:
    cmpl %edx, %esi         # if i > max_mass
    jg dp_done
    
    # Reset inner loop counter
    xorl %edi, %edi         # j = 0 (amino acid index)
    movl $.amino_count, %ebx
    
    # Get amino acid mass
    movl $.amino_masses(,%edi,4), %eax
    addl %eax, %esi         # i + mass[j]
    
    # Check if we can form this mass (i + mass[j] <= max_mass)
    cmpl %edx, %eax
    jg dp_inner_continue
    
    # dp[i] += dp[i - mass[j]]
    movl (%esi), %ecx       # temp = dp[i]
    addl (%eax), %ecx       # temp += dp[i - mass[j]]
    movl %ecx, (%esi)       # dp[i] = temp
    
dp_inner_continue:
    incl %edi               # j++
    cmpl $.amino_count, %edi
    jl dp_outer_loop
    
    jmp dp_outer_loop       # Continue with next i
    
dp_done:
    # Result is in dp[mass]
    # Assuming mass is passed as parameter or in a register
    # For demonstration, let's say we want to compute for mass = 100
    movl $100, %eax         # Example mass
    movl (%eax), %ebx       # Result = dp[mass]
    
    # Exit program
    movl $1, %eax           # sys_exit
    movl $0, %ebx           # exit status
    int $0x80

# Alternative cleaner approach for specific problem:
# Since we know the exact input format from Rosalind,
# here's a more direct implementation:

compute_peptides:
    # Input: mass in %eax
    # Output: number of peptides in %ebx
    
    # Initialize dp array (size 1001 to handle masses up to 1000)
    movl $1001, %ecx
    xorl %esi, %esi         # i = 0
    movl $1, (%esi)         # dp[0] = 1
    
    # Fill dp array from 1 to mass
    movl $1, %edi           # i = 1
    
dp_loop:
    cmpl %eax, %edi         # if i > mass
    jg dp_done
    
    # For each amino acid
    xorl %esi, %esi         # j = 0 (amino acid index)
    
inner_loop:
    # Get amino acid mass
    movl $.amino_masses(,%esi,4), %ecx
    
    # Check if we can use this amino acid
    cmpl %edi, %ecx
    jg skip_amino
    
    # dp[i] += dp[i - mass[j]]
    subl %ecx, %edi         # i - mass[j]
    movl (%edi), %edx       # dp[i - mass[j]]
    addl %edx, %eax         # dp[i] += dp[i - mass[j]]
    movl %eax, (%edi)       # Store back to dp[i]
    
skip_amino:
    incl %esi
    cmpl $20, %esi          # 20 amino acids
    jl inner_loop
    
    incl %edi               # i++
    jmp dp_loop
    
dp_done:
    movl %eax, %ebx         # Return result
    ret
```

## Key Concepts

1. **Dynamic Programming**: Build up solutions for smaller masses to solve larger ones
2. **State Definition**: `dp[i]` = number of ways to achieve mass `i`
3. **Base Case**: `dp[0] = 1` (one way to make mass 0 - empty peptide)
4. **Transition**: For each amino acid mass, add the number of ways to form the remaining mass

## Time Complexity
- **Time**: O(mass × number_of_amino_acids)
- **Space**: O(mass)

## Example Usage
For mass = 100, the program would return the number of different peptides that sum to exactly 100 daltons.

The assembly code above provides a framework for solving this problem using dynamic programming principles in x86 assembly language.