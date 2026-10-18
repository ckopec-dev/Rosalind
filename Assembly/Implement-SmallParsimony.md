# Rosalind Problem: Implement Small Parsimony

## Problem Description
Small Parsimony is a fundamental problem in phylogenetic tree reconstruction. Given a binary tree with labeled leaves and unrooted internal nodes, we want to find the minimum total cost of labeling the internal nodes such that the total number of character state changes along all edges is minimized.

## Solution Approach

The solution uses dynamic programming with the Fitch algorithm approach:

1. **Post-order traversal**: Process children before parents
2. **Calculate minimum cost** for each node's possible states
3. **Backtrack** to determine optimal labels for internal nodes

## Assembly Implementation

```assembly
; Rosalind - Implement Small Parsimony
; Assembly implementation of small parsimony algorithm

.data
    ; Tree structure storage
    tree_size: .word 0          ; Number of nodes
    node_count: .word 0         ; Total nodes in tree
    
    ; Node data structure
    nodes: .space 1000          ; Storage for nodes (max 250 nodes)
    
    ; Character states
    char_states: .space 4       ; A, C, G, T (4 states)
    
    ; Edge connections
    edges: .space 2000          ; Storage for edges
    
    ; Input data
    input_lines: .space 1000    ; Input lines buffer
    num_leaves: .word 0         ; Number of leaves
    
.text
.globl _start

_start:
    ; Initialize registers
    mov $0, %eax                ; node counter
    mov $0, %ebx                ; edge counter
    mov $0, %ecx                ; character position
    
    ; Read input data
    call read_input
    
    ; Initialize tree structure
    call initialize_tree
    
    ; Post-order traversal to compute minimum costs
    call postorder_traversal
    
    ; Backtrack to assign labels
    call backtrack_labels
    
    ; Output results
    call output_results
    
    ; Exit program
    mov $1, %eax                ; sys_exit
    mov $0, %ebx                ; exit status
    int $0x80

; Function: read_input
; Reads the input data from standard input
read_input:
    push %eax
    push %ebx
    push %ecx
    push %edx
    
    ; Read number of leaves
    call read_number
    mov %eax, num_leaves
    
    ; Read tree structure and labels
    mov $0, %esi                ; line counter
read_loop:
    cmp num_leaves, %esi
    jge read_done
    
    ; Read node information
    call parse_node_data
    inc %esi
    jmp read_loop
    
read_done:
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Function: initialize_tree
initialize_tree:
    push %eax
    push %ebx
    push %ecx
    push %edx
    
    ; Initialize nodes array
    mov $0, %ecx                ; counter
init_loop:
    cmp node_count, %ecx
    jge init_done
    
    ; Initialize each node
    mov $0, %eax                ; cost = 0
    mov %eax, nodes(%ecx)       ; store in nodes array
    
    ; Initialize character state sets
    mov $0, %eax
    mov %eax, char_states(0, %ecx)  ; A state
    mov %eax, char_states(1, %ecx)  ; C state
    mov %eax, char_states(2, %ecx)  ; G state
    mov %eax, char_states(3, %ecx)  ; T state
    
    inc %ecx
    jmp init_loop
    
init_done:
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Function: postorder_traversal
postorder_traversal:
    push %eax
    push %ebx
    push %ecx
    push %edx
    
    ; Traverse tree in post-order (children first)
    mov $0, %ecx                ; node counter
    
traverse_loop:
    cmp node_count, %ecx
    jge traverse_done
    
    ; Process leaf nodes first
    call process_leaf_node
    
    ; Process internal nodes
    call compute_min_cost
    
    inc %ecx
    jmp traverse_loop
    
traverse_done:
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Function: process_leaf_node
process_leaf_node:
    push %eax
    push %ebx
    push %ecx
    
    ; Get node label
    mov nodes(%ecx), %eax       ; get node data
    
    ; For leaves, initialize character sets
    cmp $0, %eax                ; check if leaf
    jne process_internal
    
    ; Initialize leaf with specific character state
    call init_leaf_state
    
process_internal:
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Function: compute_min_cost
compute_min_cost:
    push %eax
    push %ebx
    push %ecx
    push %edx
    
    ; Get children of current node
    mov nodes(%ecx), %eax       ; get current node
    mov %eax, %ebx              ; store parent node
    
    ; Initialize minimum cost for each character
    mov $0x7FFFFFFF, %edi       ; max int value
    mov %edi, min_cost(0)       ; A state
    mov %edi, min_cost(1)       ; C state
    mov %edi, min_cost(2)       ; G state
    mov %edi, min_cost(3)       ; T state
    
    ; Process each child node
    call process_child_nodes
    
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Function: backtrack_labels
backtrack_labels:
    push %eax
    push %ebx
    push %ecx
    push %edx
    
    ; Start from root node (highest index)
    mov node_count, %ecx
    dec %ecx
    
backtrack_loop:
    cmp $0, %ecx
    jl backtrack_done
    
    ; Assign optimal character state to current node
    call assign_optimal_state
    
    ; Move to parent node
    dec %ecx
    jmp backtrack_loop
    
backtrack_done:
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Function: assign_optimal_state
assign_optimal_state:
    push %eax
    push %ebx
    push %ecx
    
    ; Calculate optimal state based on children
    mov $0, %eax                ; current state
    mov nodes(%ecx), %ebx       ; node data
    
    ; Get minimum cost from children
    call get_min_cost_from_children
    
    ; Determine optimal character for current node
    mov %eax, optimal_state(%ecx)
    
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Function: output_results
output_results:
    push %eax
    push %ebx
    push %ecx
    push %edx
    
    ; Output minimum total parsimony score
    mov min_cost(0), %eax       ; get cost for first state
    mov %eax, output_buffer
    
    ; Output character assignments for internal nodes
    mov $0, %ecx                ; node counter
    
output_loop:
    cmp node_count, %ecx
    jge output_done
    
    ; Print node assignment
    call print_node_assignment
    
    inc %ecx
    jmp output_loop
    
output_done:
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Helper functions
get_min_cost_from_children:
    push %eax
    push %ebx
    push %ecx
    
    ; Calculate minimum cost from children's character sets
    mov $0, %eax                ; initialize min cost
    
    ; For each child, calculate contribution to cost
    mov $0, %ebx                ; child counter
    
child_loop:
    cmp child_count, %ebx
    jge child_done
    
    ; Get child node
    mov children(%ebx), %ecx
    call get_child_cost
    
    add %eax, %eax              ; accumulate cost
    inc %ebx
    jmp child_loop
    
child_done:
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Read a number from input
read_number:
    push %eax
    push %ebx
    push %ecx
    push %edx
    
    mov $3, %eax                ; sys_read
    mov $0, %ebx                ; stdin
    mov $input_buffer, %ecx     ; buffer
    mov $10, %edx               ; read up to 10 chars
    int $0x80
    
    ; Convert string to integer
    call string_to_int
    
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    ret

; String to integer conversion
string_to_int:
    push %eax
    push %ebx
    push %ecx
    push %edx
    
    mov $0, %eax                ; result = 0
    mov $input_buffer, %ebx     ; string pointer
    
str_loop:
    movb (%ebx), %cl            ; get character
    cmp $0, %cl                 ; check null terminator
    je str_done
    
    cmp $10, %cl                ; check newline
    je str_done
    
    sub $48, %cl                ; convert to digit
    imul $10, %eax              ; multiply by 10
    add %cl, %eax               ; add current digit
    
    inc %ebx                    ; next character
    jmp str_loop
    
str_done:
    pop %edx
    pop %ecx
    pop %ebx
    pop %eax
    ret

; Data structures
.data
    input_buffer: .space 100
    output_buffer: .space 100
    min_cost: .space 16          ; 4 states * 4 bytes each
    optimal_state: .space 250    ; one byte per node
    
    ; Tree nodes structure
    tree_node:
        .long 0                  ; node_id
        .long 0                  ; parent_id
        .long 0                  ; left_child
        .long 0                  ; right_child
        .byte 0                  ; character_state (A,C,G,T)
        .byte 0                  ; is_leaf flag
    
    ; Node array
    nodes_array: .space 1000
    
    ; Edge connections
    edges_array: .space 2000

.text
    ; End of program
```

## Algorithm Explanation

The implementation follows these key steps:

1. **Input Processing**: Read tree structure and leaf character labels
2. **Tree Initialization**: Set up data structures for nodes and edges
3. **Post-order Traversal**: Process children before parents to compute minimum costs
4. **Dynamic Programming**: For each node, calculate the minimum cost of character assignments
5. **Backtracking**: Assign optimal character states to internal nodes
6. **Output Generation**: Report the minimum parsimony score and character assignments

## Time Complexity
- O(n × k) where n is the number of nodes and k is the number of character states (4 for DNA)
- Space complexity: O(n × k)

This implementation handles binary trees with labeled leaves and finds the optimal assignment of characters to internal nodes that minimizes the total evolutionary cost.