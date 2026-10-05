# Find a Profile-most Probable k-mer in a String - Assembly Solution

## Problem Understanding

Given a string, a profile matrix (k × 4), and a k-mer length k, we need to find the k-mer that has the highest probability according to the profile matrix.

## Assembly Implementation

```assembly
.data
    dna_string: .ascii "ACGTACGTACGT"
    dna_length: .long 12
    k: .long 3
    profile_matrix: .float 0.2, 0.3, 0.4, 0.1   # A, C, G, T probabilities for position 1
                    .float 0.1, 0.2, 0.3, 0.4   # A, C, G, T probabilities for position 2
                    .float 0.4, 0.1, 0.2, 0.3   # A, C, G, T probabilities for position 3
    
    max_prob: .float 0.0
    best_kmer: .ascii "   "     # 3-character k-mer buffer
    temp_prob: .float 0.0

.text
.globl _start

_start:
    # Load parameters
    la $a0, dna_string      # String address
    lw $a1, dna_length      # String length
    lw $a2, k               # k-mer length
    
    # Initialize max probability to 0
    li.s $f0, 0.0           # max_prob = 0.0
    sw $f0, max_prob
    
    # Set up loop counter
    li $t0, 0               # i = 0
    add $t1, $a1, $zero     # length - k + 1 (end position)
    sub $t1, $t1, $a2      # $t1 = length - k
    addi $t1, $t1, 1       # $t1 = length - k + 1
    
loop:
    # Check if we've gone beyond the string
    bge $t0, $t1, end_loop
    
    # Extract k-mer from position i
    la $a3, dna_string     # Load string address again
    add $a3, $a3, $t0      # Point to start of k-mer
    jal extract_kmer       # Extract k-mer
    
    # Calculate probability for this k-mer
    jal calculate_probability
    
    # Compare with max probability
    l.s $f2, max_prob      # Load current max
    c.s $f1, $f2           # Compare calculated vs max
    bc1t update_max        # If calculated > max, update
    
    # Increment counter
    addi $t0, $t0, 1       # i++
    j loop
    
update_max:
    # Update max probability and best k-mer
    sw $f1, max_prob       # Store new max probability
    
    # Save the k-mer to best_kmer buffer
    la $a3, best_kmer      # Load buffer address
    li $t2, 0              # index = 0
    
copy_loop:
    beq $t2, $a2, copy_done
    lb $t3, ($a3)          # Load character from k-mer
    sb $t3, ($a4)          # Store to best_kmer buffer
    addi $a3, $a3, 1       # Move to next char in k-mer
    addi $a4, $a4, 1       # Move to next char in buffer
    addi $t2, $t2, 1       # increment index
    j copy_loop
    
copy_done:
    # Increment counter and continue loop
    addi $t0, $t0, 1
    j loop

end_loop:
    # Return best k-mer in best_kmer buffer
    la $a0, best_kmer      # Load address of result
    li $v0, 10             # Exit system call
    syscall

extract_kmer:
    # Extract k-mer from string starting at $a3 (address)
    # Return k-mer in $a4 or similar register
    addi $sp, $sp, -12     # Allocate stack space for 3 chars
    li $t0, 0              # index = 0
    
extract_loop:
    beq $t0, $a2, extract_done
    lb $t1, ($a3)          # Load character from string
    sb $t1, ($sp)          # Store to stack (k-mer)
    addi $a3, $a3, 1       # Next position in string
    addi $sp, $sp, 1       # Move stack pointer
    addi $t0, $t0, 1       # Increment index
    j extract_loop
    
extract_done:
    addi $sp, $sp, -12     # Restore stack pointer
    jr $ra                 # Return

calculate_probability:
    # Calculate probability of k-mer using profile matrix
    # Input: k-mer in $a4, profile_matrix in memory
    # Output: probability in $f1
    
    li $t0, 0              # position index = 0
    li.s $f1, 1.0          # Initialize probability to 1.0
    
calc_loop:
    beq $t0, $a2, calc_done # If all positions processed
    
    # Get character from k-mer at position $t0
    la $a3, ($a4)          # Load k-mer address
    add $a3, $a3, $t0      # Point to current position
    lb $t1, ($a3)          # Load character
    
    # Convert character to index (A=0, C=1, G=2, T=3)
    li $t2, 65             # ASCII 'A'
    sub $t3, $t1, $t2      # Get difference
    li $t4, 1              # A is at index 0, C at 1, G at 2, T at 3
    beq $t3, $zero, pos_a
    li $t5, 67             # ASCII 'C'
    sub $t6, $t1, $t5      # Get difference for C
    beq $t6, $zero, pos_c
    li $t7, 71             # ASCII 'G'  
    sub $t8, $t1, $t7      # Get difference for G
    beq $t8, $zero, pos_g
    
    # T case - index = 3
    li $t9, 84             # ASCII 'T'
    sub $t10, $t1, $t9     # Difference for T
    beq $t10, $zero, pos_t
    
pos_a:
    li $t2, 0              # Index for A
    j get_prob
    
pos_c:
    li $t2, 1              # Index for C
    j get_prob
    
pos_g:
    li $t2, 2              # Index for G
    j get_prob
    
pos_t:
    li $t2, 3              # Index for T
    
get_prob:
    # Calculate profile matrix address
    la $a3, profile_matrix # Load base address of profile
    mul $t3, $t0, 4        # Row = position * 4 (4 nucleotides)
    add $a3, $a3, $t3      # Add row offset
    mul $t4, $t2, 4        # Column = nucleotide index * 4 (float size)
    add $a3, $a3, $t4      # Add column offset
    
    l.s $f2, ($a3)         # Load probability value
    
    # Multiply with current probability
    mul.s $f1, $f1, $f2    # Update probability
    
    addi $t0, $t0, 1       # Increment position
    j calc_loop
    
calc_done:
    jr $ra                 # Return
```

## Key Algorithm Steps

1. **Extract all possible k-mers** from the DNA string
2. **Calculate probability** for each k-mer using the profile matrix
3. **Compare probabilities** and keep track of the maximum
4. **Return the best k-mer**

## Complexity Analysis

- **Time Complexity**: O(n × k) where n is string length and k is k-mer size
- **Space Complexity**: O(k) for storing current k-mer and result

## Assumptions

- Profile matrix is stored as 4 columns (A, C, G, T) by k rows
- DNA string contains only valid nucleotides (A, C, G, T)
- Input parameters are properly initialized in data section

This implementation provides a basic framework that can be adapted for specific assembly architectures and calling conventions.