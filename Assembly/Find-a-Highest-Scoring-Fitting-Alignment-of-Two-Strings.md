# Finding a Highest-Scoring Fitting Alignment of Two Strings

## Problem Understanding

The fitting alignment problem finds the highest-scoring alignment where the first string is a prefix and the second string is a suffix of the alignment, but the first string doesn't have to be aligned from position 0.

## Solution in Assembly (x86-64)

```assembly
.section .data
    # Scoring parameters
    match_score:    .long 1
    mismatch_score: .long -1
    gap_penalty:    .long -1

.section .text
    .global fitting_alignment

fitting_alignment:
    # Arguments: rdi = text1 (first string), rsi = text2 (second string)
    # Returns: highest score in rax
    
    # Save registers
    push %rbp
    push %rbx
    push %r12
    push %r13
    push %r14
    push %r15
    
    # Get string lengths
    mov $0, %rcx          # i counter for text1
    mov $0, %rdx          # j counter for text2
    
text1_len_loop:
    movb (%rdi,%rcx), %al
    testb %al, %al
    jz text1_end
    inc %rcx
    jmp text1_len_loop
    
text1_end:
    mov %rcx, %r8         # r8 = len1
    
text2_len_loop:
    movb (%rsi,%rdx), %al
    testb %al, %al
    jz text2_end
    inc %rdx
    jmp text2_len_loop
    
text2_end:
    mov %rdx, %r9         # r9 = len2
    
    # Allocate DP table: (len1+1) x (len2+1)
    # Each cell is 4 bytes (32-bit integer)
    mov %r8, %rax
    inc %rax              # +1 for row 0
    mov %r9, %rbx
    inc %rbx              # +1 for col 0
    mul %rbx              # size = (len1+1) * (len2+1)
    mov $4, %rcx          # bytes per cell
    mul %rcx              # total bytes needed
    
    # Allocate memory (simplified - in real code would use malloc)
    # For this example, we'll assume sufficient stack space
    sub %rax, %rsp        # Allocate space on stack
    
    # Initialize DP table
    # First row: 0, -1, -2, -3, ...
    mov $0, %r10          # i = 0 (first row)
    mov $0, %r11          # j = 0 (column counter)
    
init_first_row:
    cmp %r9, %r11
    jg init_first_row_done
    movl $0, (%rsp,%r11,4)  # DP[0][j] = 0
    inc %r11
    jmp init_first_row
    
init_first_row_done:
    # First column: 0, -1, -2, -3, ...
    mov $1, %r11          # j = 1 (skip first row)
    
init_first_col:
    cmp %r8, %r11
    jg init_first_col_done
    
    # DP[i][0] = i * gap_penalty
    mov %r11, %rax
    neg %rax              # -i
    imul dword ptr gap_penalty, %eax
    movl %eax, (%rsp,%r11,4)  # DP[i][0] = -i
    
    inc %r11
    jmp init_first_col
    
init_first_col_done:
    
    # Fill DP table using dynamic programming
    mov $1, %r12          # i = 1 (row counter)
    
dp_fill_loop_i:
    cmp %r8, %r12
    jg dp_fill_done
    
    mov $1, %r13          # j = 1 (column counter)
    
dp_fill_loop_j:
    cmp %r9, %r13
    jg dp_fill_next_row
    
    # Calculate DP[i][j]
    # Get characters
    movb (%rdi,%r12,1), %al     # text1[i-1]
    movb (%rsi,%r13,1), %bl     # text2[j-1]
    
    # Match/mismatch score
    cmpb %al, %bl
    je match_score_case
    
    # Mismatch case
    movl mismatch_score, %ecx
    jmp calculate_dp
    
match_score_case:
    movl match_score, %ecx
    
calculate_dp:
    # Get three possible scores:
    # 1. DP[i-1][j] + gap_penalty (deletion)
    # 2. DP[i][j-1] + gap_penalty (insertion)  
    # 3. DP[i-1][j-1] + match/mismatch score (match/mismatch)
    
    # DP[i-1][j] + gap_penalty
    mov %r12, %rax
    dec %rax              # i-1
    mov %r13, %rbx
    imul %r9, %rax        # row offset
    add %rbx, %rax        # column offset
    add %rsp, %rax        # address of DP[i-1][j]
    movl (%rax), %eax
    addl gap_penalty, %eax
    
    mov %r14, %edx        # Save this score in %edx
    
    # DP[i][j-1] + gap_penalty
    mov %r12, %rax
    imul %r9, %rax        # row offset  
    add %r13, %rax        # column offset
    dec %rax              # j-1 (adjust for 0-indexed)
    add %rsp, %rax        # address of DP[i][j-1]
    movl (%rax), %eax
    addl gap_penalty, %eax
    
    mov %r15, %ecx        # Save this score in %ecx
    
    # DP[i-1][j-1] + match/mismatch score  
    mov %r12, %rax
    dec %rax              # i-1
    imul %r9, %rax        # row offset
    add %r13, %rax        # column offset
    dec %rax              # j-1 (adjust for 0-indexed)
    add %rsp, %rax        # address of DP[i-1][j-1]
    movl (%rax), %eax
    
    addl %ecx, %eax       # Add match/mismatch score
    
    # Find maximum of three scores
    mov %edx, %ebx
    cmp %ebx, %eax
    jge max1
    mov %eax, %ebx
    
max1:
    cmp %ebx, %ecx
    jge max2
    mov %ecx, %ebx
    
max2:
    # Store result in DP[i][j]
    mov %r12, %rax
    imul %r9, %rax        # row offset
    add %r13, %rax        # column offset
    add %rsp, %rax        # address of DP[i][j]
    movl %ebx, (%rax)
    
dp_fill_next_row:
    inc %r13
    jmp dp_fill_loop_j
    
dp_fill_done:
    inc %r12
    jmp dp_fill_loop_i
    
    # Find maximum score in last row (fitting alignment constraint)
    mov $0, %r14          # max_score = 0
    
    # Check last row for maximum value
    mov %r9, %r15         # j = len2 (last column)
    
find_max_last_row:
    cmp $0, %r15
    jl find_max_done
    
    mov %r8, %rax
    imul %r9, %rax        # row offset for last row
    add %r15, %rax        # column offset
    add %rsp, %rax        # address of DP[len1][j]
    movl (%rax), %ebx
    
    cmp %r14, %ebx
    jge skip_update
    mov %ebx, %r14        # Update max score
    
skip_update:
    dec %r15
    jmp find_max_last_row
    
find_max_done:
    
    # Clean up stack
    add %rax, %rsp        # Deallocate memory
    
    # Restore registers
    pop %r15
    pop %r14
    pop %r13
    pop %r12
    pop %rbx
    pop %rbp
    
    movl %r14, %eax       # Return maximum score
    ret

# Alternative simpler approach for finding actual alignment
fitting_alignment_with_traceback:
    # Same initialization as above...
    
    # After filling DP table, do traceback from last row
    # This is more complex but provides the actual alignment
    
    # Return just the score for now (simplified)
    ret
```

## Key Concepts

1. **DP Table Construction**: 
   - Create a table of size (len1+1) × (len2+1)
   - Initialize first row with zeros (allowing gaps at start)
   - Initialize first column with gap penalties

2. **Recurrence Relation**:
   ```
   DP[i][j] = max(
       DP[i-1][j] + gap_penalty,     // deletion
       DP[i][j-1] + gap_penalty,     // insertion  
       DP[i-1][j-1] + match_score    // match
   )
   ```

3. **Fitting Constraint**:
   - The alignment must start at the beginning of text1 (first row)
   - But can end anywhere in text2 (last row)

## Time and Space Complexity

- **Time Complexity**: O(m×n) where m and n are string lengths
- **Space Complexity**: O(m×n) for the DP table

## Usage Notes

This assembly implementation:
- Assumes null-terminated strings
- Uses stack allocation for DP table (simplified)
- Returns only the maximum score
- For actual alignment, additional traceback would be needed

The fitting alignment allows for local alignment where one string is a prefix and the other is a suffix, which is useful in bioinformatics applications like sequence alignment.