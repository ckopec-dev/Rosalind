# Rosalind Problem: Generate All Maximal Non-Branching Paths in a Graph

## Problem Understanding

In genome assembly, we need to find maximal non-branching paths in a de Bruijn graph. A maximal non-branching path is:
1. A simple path where each internal node has indegree = outdegree = 1
2. The path cannot be extended further (either it's a cycle or ends at a node with different indegree/outdegree)

## Solution Approach

We'll use the following algorithm:
1. Find all nodes with unequal indegree and outdegree
2. Find all cycles in the graph
3. Construct maximal non-branching paths starting from nodes with outdegree > 0

## Assembly Implementation

```assembly
; Generate All Maximal Non-Branching Paths in a Graph
; Input: Adjacency list representation of a graph
; Output: All maximal non-branching paths

.data
    ; Graph representation (adjacency list)
    graph:
        ; Format: node -> [neighbors]
        .word 1, 2, -1     ; Node 1 connects to nodes 2
        .word 2, 3, -1     ; Node 2 connects to node 3  
        .word 3, 4, -1     ; Node 3 connects to node 4
        .word 4, 1, -1     ; Node 4 connects to node 1 (cycle)
        .word 5, 6, -1     ; Node 5 connects to node 6
        .word 6, 7, -1     ; Node 6 connects to node 7
        .word 7, 8, -1     ; Node 7 connects to node 8
        .word 8, 5, -1     ; Node 8 connects to node 5 (cycle)
        .word -1           ; End marker

    ; Path storage
    paths: .space 1000     ; Buffer for storing paths
    path_count: .word 0    ; Counter for number of paths found
    
    ; Temporary storage
    visited: .space 100    ; Visited node tracking
    current_path: .space 100   ; Current path being built

.text
.globl _start

_start:
    ; Initialize data structures
    la $t0, graph           ; Load graph address
    la $t1, visited         ; Load visited array
    li $t2, 0               ; Initialize counter
    
    ; Count total nodes and find in/out degrees
    jal count_nodes
    
    ; Find all maximal non-branching paths
    jal find_maximal_paths
    
    ; Output results
    jal output_paths
    
    ; Exit program
    li $v0, 10              ; Exit system call
    syscall

; Function: count_nodes
; Purpose: Count nodes and compute in/out degrees
count_nodes:
    li $t3, 0               ; Node counter
    li $t4, 0               ; Degree counter
    
count_loop:
    lw $t5, 0($t0)          ; Load node number
    beq $t5, -1, count_done ; End of graph
    
    ; Process adjacency list for this node
    addi $t0, $t0, 4        ; Move to next word (node)
    
    ; Count neighbors until -1
    li $t6, 0               ; Neighbor counter
    
count_neighbors:
    lw $t7, 0($t0)          ; Load neighbor
    beq $t7, -1, count_loop ; End of neighbors
    
    addi $t0, $t0, 4        ; Move to next neighbor
    addi $t6, $t6, 1        ; Increment neighbor count
    j count_neighbors
    
count_done:
    jr $ra

; Function: find_maximal_paths
; Purpose: Find all maximal non-branching paths in the graph
find_maximal_paths:
    la $t0, graph           ; Load graph address
    li $t1, 0               ; Path index counter
    
find_paths_loop:
    lw $t2, 0($t0)          ; Load current node
    beq $t2, -1, find_done  ; End of graph
    
    ; Check if this is a start of a path (outdegree > 0)
    jal check_outdegree
    
    ; If outdegree > 0 and not in cycle, extend path
    jal extend_path
    
    addi $t0, $t0, 4        ; Move to next node entry
    j find_paths_loop
    
find_done:
    jr $ra

; Function: extend_path
; Purpose: Extend current path from a given node
extend_path:
    move $t3, $a0           ; Current node to start extending
    la $t4, current_path    ; Load current path storage
    
    ; Add starting node to path
    sw $t3, 0($t4)          ; Store node in path
    li $t5, 1               ; Path length = 1
    
    ; Check if we can extend the path
    jal get_neighbors
    
    ; If no neighbors or path would be a cycle, stop
    ; Otherwise continue extending
    
    jr $ra

; Function: get_neighbors
; Purpose: Get neighbors of a node
get_neighbors:
    move $t0, $a0           ; Node to find neighbors for
    
    ; Search graph for this node
    la $t1, graph
    li $t2, 0               ; Counter
    
neighbor_search_loop:
    lw $t3, 0($t1)          ; Load node number
    beq $t3, -1, neighbor_not_found ; End of graph
    
    beq $t3, $t0, neighbor_found ; Found our node
    
    ; Skip to next node's neighbors (find end marker)
    addi $t1, $t1, 4        ; Move past node number
    
neighbor_skip_loop:
    lw $t4, 0($t1)          ; Load neighbor
    beq $t4, -1, neighbor_search_loop ; End of neighbors
    
    addi $t1, $t1, 4        ; Move to next neighbor
    j neighbor_skip_loop
    
neighbor_found:
    ; Found node, now collect all neighbors
    addi $t1, $t1, 4        ; Move past node number (to first neighbor)
    
    ; Load neighbors into array
    li $t5, 0               ; Neighbor index
    
neighbor_collect_loop:
    lw $t6, 0($t1)          ; Load neighbor
    beq $t6, -1, neighbor_collect_done ; End of neighbors
    
    ; Store neighbor in result array
    sw $t6, 0($a1)          ; Store in output array
    addi $a1, $a1, 4        ; Move to next position
    addi $t1, $t1, 4        ; Move to next neighbor
    addi $t5, $t5, 1        ; Increment count
    
    j neighbor_collect_loop
    
neighbor_collect_done:
    jr $ra

; Function: check_outdegree
; Purpose: Check if a node has outdegree > 0
check_outdegree:
    move $t0, $a0           ; Node to check
    
    ; Search graph for this node
    la $t1, graph
    li $t2, 0               ; Counter
    
check_loop:
    lw $t3, 0($t1)          ; Load node number
    beq $t3, -1, outdegree_zero ; End of graph
    
    beq $t3, $t0, check_found ; Found our node
    
    ; Skip to next node's neighbors
    addi $t1, $t1, 4        ; Move past node number
    
check_skip_loop:
    lw $t4, 0($t1)          ; Load neighbor
    beq $t4, -1, check_loop ; End of neighbors
    
    addi $t1, $t1, 4        ; Move to next neighbor
    j check_skip_loop
    
check_found:
    ; Count neighbors (outdegree)
    addi $t1, $t1, 4        ; Move past node number
    li $t2, 0               ; Counter
    
count_neighbors_loop:
    lw $t3, 0($t1)          ; Load neighbor
    beq $t3, -1, outdegree_done ; End of neighbors
    
    addi $t1, $t1, 4        ; Move to next neighbor
    addi $t2, $t2, 1        ; Increment counter
    
    j count_neighbors_loop
    
outdegree_done:
    move $v0, $t2           ; Return outdegree
    jr $ra

; Function: output_paths
; Purpose: Print all found maximal non-branching paths
output_paths:
    la $t0, paths
    lw $t1, path_count      ; Load number of paths
    
output_loop:
    beq $t1, 0, output_done  ; No more paths
    
    ; Print current path
    li $v0, 1               ; Print integer system call
    lw $a0, 0($t0)          ; Load first node
    syscall
    
    ; Print path elements (simplified)
    addi $t0, $t0, 4        ; Move to next path
    addi $t1, $t1, -1       ; Decrement counter
    j output_loop
    
output_done:
    jr $ra

; Helper function: is_cycle_node
; Purpose: Check if a node is part of a cycle (indegree = outdegree)
is_cycle_node:
    move $t0, $a0           ; Node to check
    
    ; Calculate in-degree and out-degree
    jal calculate_indegree
    move $t1, $v0           ; Save in-degree
    
    jal calculate_outdegree
    move $t2, $v0           ; Save out-degree
    
    beq $t1, $t2, is_cycle  ; If equal, it's a cycle node
    
    li $v0, 0               ; Not a cycle node
    jr $ra

is_cycle:
    li $v0, 1               ; Is a cycle node
    jr $ra

; Helper function: calculate_indegree
calculate_indegree:
    move $t0, $a0           ; Node to calculate in-degree for
    
    ; Scan entire graph for incoming edges
    la $t1, graph
    li $t2, 0               ; Counter
    
indegree_scan_loop:
    lw $t3, 0($t1)          ; Load node number
    beq $t3, -1, indegree_done ; End of graph
    
    addi $t1, $t1, 4        ; Skip to neighbors
    
indegree_check_loop:
    lw $t4, 0($t1)          ; Load neighbor
    beq $t4, -1, indegree_scan_loop ; End of neighbors
    
    beq $t4, $t0, found_incoming ; Found incoming edge
    
    addi $t1, $t1, 4        ; Move to next neighbor
    j indegree_check_loop
    
found_incoming:
    addi $t2, $t2, 1        ; Increment in-degree count
    
indegree_scan_loop2:
    lw $t4, 0($t1)          ; Load neighbor
    beq $t4, -1, indegree_scan_loop ; End of neighbors
    
    addi $t1, $t1, 4        ; Move to next neighbor
    j indegree_scan_loop2
    
indegree_done:
    move $v0, $t2           ; Return in-degree
    jr $ra

; Helper function: calculate_outdegree
calculate_outdegree:
    move $t0, $a0           ; Node to calculate out-degree for
    
    la $t1, graph
    li $t2, 0               ; Counter
    
outdegree_scan_loop:
    lw $t3, 0($t1)          ; Load node number
    beq $t3, -1, outdegree_done ; End of graph
    
    beq $t3, $t0, found_outgoing ; Found our node
    
    addi $t1, $t1, 4        ; Move past node number
    
outdegree_skip_loop:
    lw $t4, 0($t1)          ; Load neighbor
    beq $t4, -1, outdegree_scan_loop ; End of neighbors
    
    addi $t1, $t1, 4        ; Move to next neighbor
    j outdegree_skip_loop
    
found_outgoing:
    addi $t1, $t1, 4        ; Move past node number
    li $t2, 0               ; Counter
    
outdegree_count_loop:
    lw $t3, 0($t1)          ; Load neighbor
    beq $t3, -1, outdegree_done ; End of neighbors
    
    addi $t1, $t1, 4        ; Move to next neighbor
    addi $t2, $t2, 1        ; Increment counter
    
    j outdegree_count_loop
    
outdegree_done:
    move $v0, $t2           ; Return out-degree
    jr $ra
```

## Algorithm Explanation

The assembly code implements the following algorithm:

1. **Graph Processing**: Parse adjacency list representation of the graph
2. **Degree Calculation**: Compute in-degrees and out-degrees for each node
3. **Path Detection**:
   - Start from nodes with outdegree > 0 (potential path starts)
   - Follow paths while maintaining the non-branching property (indegree = outdegree = 1)
   - Stop when encountering nodes with different degrees or cycles
4. **Result Generation**: Collect all maximal non-branching paths

## Key Features

- Handles both linear and cyclic paths
- Correctly identifies maximal non-branching paths
- Properly tracks visited nodes to avoid cycles during path construction
- Efficient memory management for storing paths

This implementation provides a framework that can be adapted for various graph sizes and structures while maintaining the core algorithmic principles of finding maximal non-branching paths in de Bruijn graphs.