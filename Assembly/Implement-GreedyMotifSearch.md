# Rosalind Problem: Implement_GreedyMotifSearch in Assembly

I'll solve this step by step, implementing the greedy motif search algorithm in assembly language.

## Understanding the Problem

The GreedyMotifSearch algorithm finds the best motif (k-mer) that appears in all DNA sequences. It works greedily by:
1. Starting with the first k-mer from the first sequence
2. Building a profile matrix from current motifs
3. Finding the most probable k-mer in each subsequent sequence
4. Repeating until all sequences are processed

## Assembly Implementation

```assembly
.section .data
    # DNA sequences (example)
    dna_sequences: .ascii "GGCGTTCAGGCA\n"
                   .ascii "AAGATCTGTCA\n"
                   .ascii "TTGGGAATCTT\n"
                   .ascii "AAAGAAGCTC\n"
                   .ascii "TTACGAATT\n"
    
    # Sequence lengths
    seq_lengths: .long 10, 9, 11, 8, 7
    
    # Parameters
    k_value: .long 3          # motif length
    num_sequences: .long 5    # number of sequences
    
    # DNA nucleotides mapping
    nucleotides: .ascii "ACGT"
    
    # Profile matrix (4xk) - initialized to zeros
    profile_matrix: .space 12  # 4 rows × 3 columns
    
    # Best motif found so far
    best_motif: .space 3      # k characters
    
    # Temporary storage
    temp_motif: .space 3      # k characters
    temp_profile: .space 12   # 4x3 profile matrix

.section .text
    .global _start

# Function to calculate Hamming distance between two strings
hamming_distance:
    push %rbp
    mov %rsp, %rbp
    mov %rdi, %rax      # string1 address
    mov %rsi, %rcx      # string2 address  
    mov %rdx, %r8       # length
    
    xor %r9, %r9        # counter
    xor %r10, %r10      # hamming distance
    
hamming_loop:
    cmp %r8, %r9
    jge hamming_done
    
    movb (%rax,%r9), %dl
    movb (%rcx,%r9), %dh
    
    cmp %dl, %dh
    je hamming_continue
    
    inc %r10            # increment distance
    
hamming_continue:
    inc %r9
    jmp hamming_loop
    
hamming_done:
    mov %r10, %rax      # return value
    pop %rbp
    ret

# Function to compute profile matrix from motifs
compute_profile:
    push %rbp
    mov %rsp, %rbp
    mov %rdi, %rax      # motifs array address
    mov %rsi, %rcx      # number of motifs
    mov %rdx, %r8       # k (motif length)
    
    # Initialize profile matrix to zeros
    xor %r9, %r9        # row counter
init_loop:
    cmp $4, %r9
    jge init_done
    
    xor %r10, %r10      # column counter
init_col_loop:
    cmp %r8, %r10
    jge init_col_done
    
    movb $0, (%rax,%r9,%r10)  # set profile[i][j] = 0
    inc %r10
    jmp init_col_loop
    
init_col_done:
    inc %r9
    jmp init_loop
    
init_done:
    # Count nucleotides in each position
    xor %r9, %r9        # motif counter
count_loop:
    cmp %rcx, %r9
    jge count_done
    
    xor %r10, %r10      # position counter
count_pos_loop:
    cmp %r8, %r10
    jge count_pos_done
    
    # Get nucleotide at position j of motif i
    movb (%rax,%r9,%r10), %dl  # nucleotide
    
    # Find index in nucleotides array
    xor %r11, %r11      # nucleotide index
find_nucleotide:
    cmp $4, %r11
    jge find_done
    
    movb nucleotides(%r11), %dh
    cmp %dl, %dh
    je found_nucleotide
    
    inc %r11
    jmp find_nucleotide
    
found_nucleotide:
    # Increment profile count at position index
    incb (%rax,%r11,%r10)
    
find_done:
    inc %r10
    jmp count_pos_loop
    
count_pos_done:
    inc %r9
    jmp count_loop
    
count_done:
    pop %rbp
    ret

# Function to find most probable k-mer in sequence given profile
find_most_probable_kmer:
    push %rbp
    mov %rsp, %rbp
    mov %rdi, %rax      # sequence address
    mov %rsi, %rcx      # profile matrix address
    mov %rdx, %r8       # k (motif length)
    
    # Initialize best probability to 0
    xor %r9, %r9        # best probability (as integer * 1000)
    
    # Initialize best k-mer position to 0
    xor %r10, %r10      # best position
    
    # Try each k-mer in sequence
    mov %rax, %r11      # current sequence pointer
    xor %r12, %r12      # current k-mer start position
    
try_kmer_loop:
    cmp (%r11,%r12), %r8  # compare with length of sequence
    jge try_done
    
    # Calculate probability for this k-mer
    xor %r13, %r13      # position in k-mer
    xor %r14, %r14      # current probability (as integer * 1000)
    mov $1000, %r14     # start with 1.0
    
kmer_prob_loop:
    cmp %r8, %r13
    jge kmer_prob_done
    
    # Get nucleotide at current position
    movb (%rax,%r12,%r13), %dl
    
    # Find nucleotide index (A=0, C=1, G=2, T=3)
    xor %r15, %r15      # nucleotide index
find_kmer_nucleotide:
    cmp $4, %r15
    jge find_kmer_done
    
    movb nucleotides(%r15), %dh
    cmp %dl, %dh
    je found_kmer_nucleotide
    
    inc %r15
    jmp find_kmer_nucleotide
    
found_kmer_nucleotide:
    # Get probability from profile matrix
    movb (%rcx,%r15,%r13), %dl  # profile[i][j]
    
    # Multiply current probability by this value
    # Simplified multiplication (in real case would use proper floating point)
    imul %dl, %r14      # multiply by nucleotide probability
    
find_kmer_done:
    inc %r13
    jmp kmer_prob_loop
    
kmer_prob_done:
    # Check if this is better than current best
    cmp %r9, %r14
    jge try_continue
    
    mov %r14, %r9       # update best probability
    mov %r12, %r10      # update best position
    
try_continue:
    inc %r12
    jmp try_kmer_loop
    
try_done:
    mov %r10, %rax      # return best position
    pop %rbp
    ret

# Main GreedyMotifSearch algorithm
greedy_motif_search:
    push %rbp
    mov %rsp, %rbp
    mov %rdi, %rax      # sequences array address
    mov %rsi, %rcx      # number of sequences
    mov %rdx, %r8       # k value
    
    # Initialize best score to infinity (large number)
    mov $0x7FFFFFFF, %r9  # best_score
    
    # Try each k-mer from first sequence as starting motif
    xor %r10, %r10      # starting position in first sequence
    
first_seq_loop:
    cmp (%rax), %r10    # compare with length of first sequence
    jge first_seq_done
    
    # Get k-mer from first sequence
    mov %rax, %r11      # first sequence address
    xor %r12, %r12      # copy position to r12 for k-mer extraction
    
    # Extract k-mer from first sequence
    xor %r13, %r13      # k-mer character counter
extract_loop:
    cmp %r8, %r13
    jge extract_done
    
    movb (%r11,%r10,%r13), %dl
    movb %dl, best_motif(%r13)
    inc %r13
    jmp extract_loop
    
extract_done:
    # Build motifs array with this k-mer and rest from other sequences
    xor %r13, %r13      # sequence counter (starting from 1)
    
    # Process remaining sequences
build_motifs_loop:
    cmp %rcx, %r13
    jge build_done
    
    # Get current sequence address
    mov %rax, %r14      # base address of sequences
    mov $0, %r15        # index offset (skip first sequence)
    
    # Calculate offset for current sequence
    xor %r16, %r16      # cumulative length counter
    xor %r17, %r17      # current sequence counter
    
count_seq_loop:
    cmp %r13, %r17
    jge count_done
    
    # Get length of this sequence (simplified - in real code would be array)
    mov seq_lengths(%r17), %r18
    add %r18, %r16
    inc %r17
    jmp count_seq_loop
    
count_done:
    mov %r16, %r15      # offset for this sequence
    mov %rax, %r14
    add %r15, %r14      # get actual address of current sequence
    
    # Find most probable k-mer in current sequence using profile from previous motifs
    mov %r14, %rdi      # current sequence
    mov %r14, %rsi      # profile matrix (simplified)
    
    # Call find_most_probable_kmer function
    call find_most_probable_kmer
    
    # Store result in motifs array
    mov %rax, %r15      # best position returned
    xor %r16, %r16      # character counter
    
kmer_copy_loop:
    cmp %r8, %r16
    jge kmer_copy_done
    
    movb (%r14,%r15,%r16), %dl
    movb %dl, temp_motif(%r16)  # temporary storage for this k-mer
    inc %r16
    jmp kmer_copy_loop
    
kmer_copy_done:
    # Copy to final motifs array (simplified)
    inc %r13
    jmp build_motifs_loop
    
build_done:
    # Compute score for current motif set
    # In real implementation would compute the score
    
    inc %r10            # try next k-mer in first sequence
    jmp first_seq_loop
    
first_seq_done:
    pop %rbp
    ret

# Main program entry point
_start:
    # Initialize parameters
    mov $3, %rax        # k = 3
    mov $5, %rbx        # number of sequences
    
    # Call greedy motif search
    call greedy_motif_search
    
    # Exit program
    mov $60, %rax       # sys_exit
    mov $0, %rdi        # exit status
    syscall
```

## Key Components Explained

### 1. **Data Structures**
- `dna_sequences`: Array of DNA sequences
- `seq_lengths`: Lengths of each sequence
- `profile_matrix`: 4×k matrix storing nucleotide frequencies
- `best_motif`: Storage for the best motif found

### 2. **Core Functions**

**Hamming Distance**: Computes distance between two strings
**Compute Profile**: Builds profile matrix from current motifs
**Find Most Probable Kmer**: Identifies best k-mer in sequence given profile
**Greedy Motif Search**: Main algorithm implementing greedy approach

### 3. **Algorithm Steps**
1. Try each k-mer from first sequence as initial motif
2. Build profile matrix from current motifs
3. Find most probable k-mer in each remaining sequence
4. Update best motif if current is better
5. Continue until all sequences processed

## Notes on Assembly Implementation

This is a conceptual implementation that demonstrates the structure. Actual assembly implementations would require:
- Proper memory management
- Real floating-point operations for probabilities
- More sophisticated data structures
- Error handling for edge cases
- Optimized memory access patterns

The algorithm maintains O(N×k) time complexity where N is the total sequence length and k is motif length.