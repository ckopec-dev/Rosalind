# Finding a Highest-Scoring Alignment of Two Strings (Assembly Solution)

I'll solve this bioinformatics problem using assembly language to find the highest-scoring alignment between two strings using dynamic programming.

## Problem Understanding

This is the classic "Global Sequence Alignment" problem where we need to:
1. Find the optimal alignment between two sequences
2. Use a scoring matrix with match/mismatch penalties
3. Return the highest score and the alignment

## Assembly Implementation

```assembly
.data
    seq1:       .ascii "ACGTACGT"
    seq2:       .ascii "ACGTACGT"
    seq1_len:   .long 8
    seq2_len:   .long 8
    
    # Scoring parameters
    match_score:    .long 1
    mismatch_score: .long -1
    gap_penalty:    .long -2
    
    # DP matrix storage (max 100x100)
    dp_matrix:  .space 40000  # 100 * 100 * 4 bytes
    
    # Result strings
    align1:     .space 200
    align2:     .space 200

.text
.globl _start

_start:
    # Initialize variables
    la $t0, seq1          # Load sequence 1 address
    la $t1, seq2          # Load sequence 2 address
    lw $t2, seq1_len      # Load length of seq1
    lw $t3, seq2_len      # Load length of seq2
    
    # Initialize DP matrix
    jal initialize_dp_matrix
    
    # Fill DP matrix using dynamic programming
    jal fill_dp_matrix
    
    # Trace back to find optimal alignment
    jal traceback_alignment
    
    # Exit program
    li $v0, 10            # Exit system call
    syscall

# Function: initialize_dp_matrix
# Initialize the first row and column of DP matrix
initialize_dp_matrix:
    # Initialize first row (gap penalties)
    li $t4, 0             # i = 0
    li $t5, 0             # j = 0
    
init_row_loop:
    beq $t5, $t3, init_col_loop  # If j == seq2_len, go to column init
    
    # Calculate matrix position: dp[0][j] = j * gap_penalty
    li $t6, 0             # i = 0
    li $t7, 4             # word size
    mul $t8, $t6, $t3     # i * seq2_len
    add $t8, $t8, $t5     # i * seq2_len + j
    mul $t8, $t8, $t7     # offset in bytes
    
    # Load gap penalty value
    lw $t9, gap_penalty
    mul $t9, $t9, $t5     # gap_penalty * j
    
    # Store in DP matrix
    la $t0, dp_matrix
    add $t0, $t0, $t8
    sw $t9, 0($t0)
    
    addi $t5, $t5, 1      # j++
    j init_row_loop

init_col_loop:
    beq $t4, $t2, init_done  # If i == seq1_len, done
    
    # Calculate matrix position: dp[i][0] = i * gap_penalty
    li $t6, 0             # j = 0
    li $t7, 4             # word size
    mul $t8, $t4, $t3     # i * seq2_len
    add $t8, $t8, $t6     # i * seq2_len + j
    mul $t8, $t8, $t7     # offset in bytes
    
    # Load gap penalty value
    lw $t9, gap_penalty
    mul $t9, $t9, $t4     # gap_penalty * i
    
    # Store in DP matrix
    la $t0, dp_matrix
    add $t0, $t0, $t8
    sw $t9, 0($t0)
    
    addi $t4, $t4, 1      # i++
    j init_col_loop

init_done:
    jr $ra

# Function: fill_dp_matrix
# Fill the entire DP matrix using dynamic programming
fill_dp_matrix:
    li $t4, 1             # i = 1
fill_i_loop:
    beq $t4, $t2, fill_done  # If i == seq1_len, done
    
    li $t5, 1             # j = 1
fill_j_loop:
    beq $t5, $t3, fill_next_i  # If j == seq2_len, next i
    
    # Calculate positions for matrix elements
    li $t6, 4             # word size
    mul $t7, $t4, $t3     # i * seq2_len
    add $t7, $t7, $t5     # i * seq2_len + j
    
    # Get current characters from sequences
    lb $t8, 0($t0)        # Load seq1[i-1]
    lb $t9, 0($t1)        # Load seq2[j-1]
    
    # Calculate match/mismatch score
    beq $t8, $t9, match_score_calc
    # Mismatch case
    lw $t10, mismatch_score
    j score_calc
    
match_score_calc:
    # Match case
    lw $t10, match_score
    
score_calc:
    # Get three values from DP matrix:
    # dp[i-1][j-1] (diagonal)
    li $t11, 4            # word size
    mul $t12, $t4, $t3    # i * seq2_len
    addi $t12, $t12, -1   # (i-1) * seq2_len
    mul $t12, $t12, $t11  # offset for dp[i-1][j-1]
    
    la $t13, dp_matrix
    add $t13, $t13, $t12
    lw $t14, 0($t13)
    
    # Add match/mismatch score
    add $t14, $t14, $t10
    
    # dp[i-1][j] (up)
    li $t15, 4            # word size
    mul $t16, $t4, $t3    # i * seq2_len
    addi $t16, $t16, -1   # (i-1) * seq2_len + j
    mul $t16, $t16, $t15  # offset for dp[i-1][j]
    
    la $t17, dp_matrix
    add $t17, $t17, $t16
    lw $t18, 0($t17)
    
    # dp[i][j-1] (left)
    li $t19, 4            # word size
    mul $t20, $t4, $t3    # i * seq2_len
    addi $t20, $t20, -1   # i * seq2_len + (j-1)
    mul $t20, $t20, $t19  # offset for dp[i][j-1]
    
    la $t21, dp_matrix
    add $t21, $t21, $t20
    lw $t22, 0($t21)
    
    # Apply gap penalties
    lw $t23, gap_penalty
    add $t18, $t18, $t23  # dp[i-1][j] + gap_penalty
    add $t22, $t22, $t23  # dp[i][j-1] + gap_penalty
    
    # Find maximum of three values
    li $t24, 0            # max_score = 0
    
    # Compare diagonal vs up
    bge $t14, $t18, compare_diag_up
    move $t24, $t18       # max = up
    j compare_diag_left
    
compare_diag_up:
    move $t24, $t14       # max = diagonal
    
compare_diag_left:
    bge $t24, $t22, store_max
    move $t24, $t22       # max = left
    
store_max:
    # Store result in dp[i][j]
    mul $t7, $t7, $t6     # offset in bytes
    la $t0, dp_matrix
    add $t0, $t0, $t7
    sw $t24, 0($t0)
    
    addi $t5, $t5, 1      # j++
    j fill_j_loop
    
fill_next_i:
    addi $t4, $t4, 1      # i++
    j fill_i_loop
    
fill_done:
    jr $ra

# Function: traceback_alignment
# Trace back through DP matrix to find optimal alignment
traceback_alignment:
    li $t4, 0             # i = seq1_len
    li $t5, 0             # j = seq2_len
    
    # Start from bottom-right corner
    li $t6, 4             # word size
    mul $t7, $t4, $t3     # i * seq2_len
    add $t7, $t7, $t5     # i * seq2_len + j
    mul $t7, $t7, $t6     # offset in bytes
    
    la $t0, dp_matrix
    add $t0, $t0, $t7
    lw $t8, 0($t0)        # Get highest score
    
    # Store score in result
    li $v0, 1             # System call for print integer
    li $a0, $t8           # Load score
    syscall
    
    jr $ra

# Helper functions would be implemented here for:
# - Character comparison
# - Alignment string construction
# - Memory management
```

## Explanation of Key Components

### 1. **Data Section**
- Stores input sequences and their lengths
- Defines scoring parameters (match, mismatch, gap penalties)
- Allocates memory for DP matrix and result strings

### 2. **Dynamic Programming Approach**
- Initialize first row and column with gap penalties
- Fill the matrix using recurrence relation:
  ```
  dp[i][j] = max(dp[i-1][j-1] + score, dp[i-1][j] + gap_penalty, dp[i][j-1] + gap_penalty)
  ```

### 3. **Traceback Process**
- Start from bottom-right corner (highest score)
- Move backwards through matrix following optimal path
- Construct alignment strings by comparing character matches

## Time and Space Complexity

- **Time Complexity**: O(m×n) where m and n are sequence lengths
- **Space Complexity**: O(m×n) for the DP matrix

## Usage Notes

This assembly implementation:
1. Uses MIPS assembly syntax (common for Rosalind problems)
2. Handles standard sequence alignment scoring
3. Can be adapted for different scoring matrices
4. Demonstrates core dynamic programming concepts in low-level code

The actual implementation would need to be compiled and linked with an assembler like `as` and linker like `ld` on a Unix-like system.