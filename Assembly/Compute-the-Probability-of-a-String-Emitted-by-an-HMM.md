# Rosalind Problem: Compute the Probability of a String Emitted by an HMM

## Problem Understanding

We need to compute the probability that a Hidden Markov Model (HMM) emits a specific string. This is typically solved using the forward algorithm.

## Solution in Assembly Language

```assembly
; Compute the probability of a string emitted by an HMM
; Using the Forward Algorithm

.section .data
    ; HMM parameters would be loaded here
    ; State transitions, emission probabilities, initial probabilities
    
    ; Example: 3 states (A, B, C)
    ; Emission matrix (3x4) for nucleotides A,C,G,T
    ; Transition matrix (3x3) 
    ; Initial probabilities (3)

.section .text
    .global _start

_compute_probability:
    ; Input: string to be emitted, HMM parameters
    ; Output: probability of string
    
    ; Initialize forward variables
    ; For each state, compute forward probability at each position
    
    movl $string_length, %ecx        ; Loop counter (length of emitted string)
    movl $num_states, %edx           ; Number of states
    
    ; Initialize forward matrix F[0][i] = π[i] * b_i(string[0])
    call initialize_forward
    
    ; Forward algorithm loop
forward_loop:
    testl %ecx, %ecx                 ; Check if string is exhausted
    jz probability_computed
    
    ; For each state j at current position t
    movl $0, %esi                    ; State index i = 0
state_loop:
    ; Compute F[t][j] = sum over all states i of [F[t-1][i] * a_ij * b_j(string[t])]
    
    ; F[t][j] = sum_i(F[t-1][i] * a_ij * b_j(string[t]))
    
    call compute_forward_value
    
    incl %esi                        ; Next state
    cmpl %edx, %esi                  ; Check if all states processed
    jl state_loop
    
    decl %ecx                        ; Next character in string
    jmp forward_loop

compute_forward_value:
    ; Compute F[t][j] = sum_i(F[t-1][i] * a_ij * b_j(string[t]))
    
    ; Load F[t-1][i]
    movl forward_matrix(%esi), %eax  ; F[t-1][i]
    
    ; Multiply by transition probability a_ij
    movl transition_matrix(%esi, %edi, 4), %ebx  ; a_ij
    imull %ebx, %eax
    
    ; Multiply by emission probability b_j(string[t])
    movl emission_matrix(%esi, %edi, 4), %ebx   ; b_j(string[t])
    imull %ebx, %eax
    
    ; Add to running sum
    addl %eax, %ebp                  ; Accumulate in sum register
    
    ret

initialize_forward:
    ; Initialize F[0][i] = π[i] * b_i(string[0])
    movl $0, %esi                    ; State index i = 0
    
init_loop:
    ; F[0][i] = π[i] * b_i(string[0])
    
    ; Load initial probability π[i]
    movl initial_prob(%esi), %eax
    
    ; Load emission probability b_i(string[0])
    movl emission_matrix(%esi, %edi, 4), %ebx
    
    ; Multiply and store
    imull %ebx, %eax
    movl %eax, forward_matrix(%esi)
    
    incl %esi                        ; Next state
    cmpl %edx, %esi                  ; Check if all states processed
    jl init_loop
    
    ret

probability_computed:
    ; Return final probability in eax
    movl forward_matrix, %eax        ; Final result
    
    ret

; Example usage:
; Call with string pointer and HMM parameters loaded
; Function returns probability as integer (scaled)
```

## Alternative Implementation with Better Structure

```assembly
; More structured implementation of the Forward Algorithm

.section .data
    # HMM parameters would be defined here
    # States: A, B, C (3 states)
    # Alphabet: A, C, G, T (4 symbols)
    
    num_states:     .long 3
    num_symbols:    .long 4
    string_length:  .long 0
    
    # Initial probabilities π[i]
    initial_probs:  .long 0.5, 0.3, 0.2   # Example values
    
    # Transition matrix a_ij
    transition_matrix:
        .long 0.7, 0.2, 0.1   # A -> A, B, C
        .long 0.3, 0.5, 0.2   # B -> A, B, C  
        .long 0.4, 0.1, 0.5   # C -> A, B, C
    
    # Emission matrix b_i(symbol)
    emission_matrix:
        .long 0.3, 0.2, 0.3, 0.2   # State A emissions for A,C,G,T
        .long 0.1, 0.4, 0.3, 0.2   # State B emissions for A,C,G,T
        .long 0.2, 0.3, 0.2, 0.3   # State C emissions for A,C,G,T

.section .text
    .global forward_algorithm

forward_algorithm:
    ; Input: string pointer, length
    ; Output: probability
    
    pushl %ebp
    movl %esp, %ebp
    
    ; Get parameters from stack
    movl 8(%ebp), %esi     ; string pointer
    movl 12(%ebp), %ecx    ; string length
    
    ; Initialize forward matrix
    call initialize_forward_matrix
    
    ; Main forward loop
    movl $0, %edi          ; position counter = 0
    
forward_main_loop:
    cmpl %ecx, %edi        ; check if we've processed all characters
    jge end_forward
    
    ; Process current character (character at position %edi)
    call process_current_position
    
    incl %edi              ; next position
    jmp forward_main_loop

end_forward:
    ; Sum up final probabilities
    call sum_final_probabilities
    movl %eax, %esi        ; return probability
    
    popl %ebp
    ret

initialize_forward_matrix:
    ; Initialize F[0][i] = π[i] * b_i(string[0])
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
    
    movl num_states, %ecx
    movl $0, %edx          ; state index
    
init_loop:
    ; Get initial probability π[i]
    movl initial_probs(,%edx,4), %eax
    
    ; Get emission probability b_i(string[0])
    movl (%esi), %ebx      ; first character of string
    movl emission_matrix(,%edx,4), %ebx
    
    ; Multiply: F[0][i] = π[i] * b_i(char_0)
    imull %ebx, %eax
    movl %eax, forward_matrix(,%edx,4)  ; Store in forward matrix
    
    incl %edx
    decl %ecx
    jnz init_loop
    
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
    ret

process_current_position:
    ; Process position t, update F[t][j] for all states j
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
    
    movl num_states, %ecx      ; number of states
    movl $0, %edx              ; state index j = 0
    
process_state_loop:
    ; Compute F[t][j] = sum_i(F[t-1][i] * a_ij * b_j(x_t))
    
    ; Initialize sum to 0
    movl $0, %ebx
    
    ; For each previous state i
    movl num_states, %esi      ; number of states
    movl $0, %eax              ; previous state index i = 0
    
compute_sum_loop:
    ; Get F[t-1][i]
    movl forward_matrix(,%eax,4), %edi
    
    ; Multiply by transition probability a_ij
    movl transition_matrix(,%eax,4), %ebp
    imull %ebp, %edi
    
    ; Multiply by emission probability b_j(x_t)
    movl emission_matrix(,%edx,4), %ebp
    imull %ebp, %edi
    
    ; Add to sum
    addl %edi, %ebx
    
    incl %eax
    decl %esi
    jnz compute_sum_loop
    
    ; Store result F[t][j] = sum
    movl %ebx, forward_matrix(,%edx,4)
    
    incl %edx
    decl %ecx
    jnz process_state_loop
    
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
    ret

sum_final_probabilities:
    ; Sum up all final probabilities to get total probability
    pushl %eax
    pushl %ebx
    pushl %ecx
    
    movl num_states, %ecx
    movl $0, %ebx           ; sum = 0
    
final_sum_loop:
    addl forward_matrix(,%ecx,4), %ebx
    decl %ecx
    jnz final_sum_loop
    
    movl %ebx, %eax         ; return result
    
    popl %ecx
    popl %ebx
    popl %eax
    ret

; Helper function to get character from string
get_char:
    ; Input: string pointer, position
    ; Output: character value
    pushl %ebp
    movl %esp, %ebp
    
    movl 8(%ebp), %esi     ; string pointer
    movl 12(%ebp), %edi    ; position
    
    movb (%esi,%edi,1), %al ; get character at position
    
    popl %ebp
    ret
```

## Key Concepts Implemented

1. **Forward Algorithm**: Dynamic programming approach to compute emission probabilities
2. **Matrix Operations**: Handling transition and emission matrices
3. **State Tracking**: Maintaining forward variables for each state at each time step
4. **Probability Accumulation**: Multiplying probabilities along the HMM paths

## Time Complexity
- O(T × N²) where T is string length and N is number of states
- Space complexity: O(N × T) for forward matrix

This assembly implementation provides a foundation for solving the HMM probability computation problem, though actual implementation would require careful handling of floating-point arithmetic and memory management.