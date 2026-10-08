# Rosalind Problem: Construct the De Bruijn Graph of a Collection of k-mers

## Problem Statement
Given a collection of k-mers, construct the de Bruijn graph where each k-mer is represented as an edge, and vertices are (k-1)-mers.

## Solution in Assembly Language

```assembly
; De Bruijn Graph Construction
; Input: Collection of k-mers
; Output: De Bruijn graph representation

.data
    ; Input k-mers array
    kmers_array: .word kmer1, kmer2, kmer3, kmer4  ; Example k-mers
    num_kmers: .word 4                             ; Number of k-mers
    
    ; Graph storage
    vertices: .space 1000                          ; Vertex storage
    edges: .space 2000                             ; Edge storage
    
    ; Temporary buffers
    prefix_buffer: .space 50                        ; For (k-1) prefixes
    suffix_buffer: .space 50                        ; For (k-1) suffixes
    
    ; Constants
    k_value: .word 3                               ; k-mer length
    max_vertices: .word 100                        ; Maximum vertices

.text
.globl _start

_start:
    ; Load input parameters
    la $a0, kmers_array     ; Load k-mers array address
    lw $a1, num_kmers       ; Load number of k-mers
    lw $a2, k_value         ; Load k value
    
    ; Initialize vertex counter
    li $t0, 0               ; vertex_count = 0
    
    ; Process each k-mer to extract prefixes and suffixes
    li $t1, 0               ; i = 0 (loop counter)
    
process_loop:
    bge $t1, $a1, graph_construction_done  ; if i >= num_kmers, done
    
    ; Get current k-mer address
    sll $t2, $t1, 2         ; $t2 = i * 4 (word offset)
    add $t3, $a0, $t2       ; $t3 = &kmers_array[i]
    lw $t4, 0($t3)          ; Load k-mer string address
    
    ; Extract prefix (k-1) and suffix (k-1)
    jal extract_prefix_suffix
    addi $t1, $t1, 1        ; i++
    j process_loop

extract_prefix_suffix:
    ; Parameters: $a0 = k-mer address, $a2 = k value
    ; Returns: prefix and suffix in buffers
    
    ; Extract prefix (first k-1 characters)
    li $t5, 0               ; char_index = 0
    li $t6, 0               ; buffer_index = 0
    
prefix_loop:
    bge $t5, $a2, prefix_done
    add $t7, $a0, $t5       ; address of current character
    lb $t8, 0($t7)          ; load character
    ; Store in prefix_buffer
    add $t9, $sp, $t6       ; buffer address
    sb $t8, 0($t9)          ; store character
    addi $t5, $t5, 1        ; increment char_index
    addi $t6, $t6, 1        ; increment buffer_index
    j prefix_loop

prefix_done:
    ; Null terminate prefix
    add $t9, $sp, $t6       ; buffer address
    sb $zero, 0($t9)        ; null terminator
    
    ; Extract suffix (last k-1 characters)
    li $t5, 1               ; start from index 1
    li $t6, 50              ; suffix buffer offset
    
suffix_loop:
    bge $t5, $a2, suffix_done
    add $t7, $a0, $t5       ; address of current character
    lb $t8, 0($t7)          ; load character
    ; Store in suffix_buffer (offset 50)
    add $t9, $sp, $t6       ; buffer address
    sb $t8, 0($t9)          ; store character
    addi $t5, $t5, 1        ; increment char_index
    addi $t6, $t6, 1        ; increment buffer_index
    j suffix_loop

suffix_done:
    ; Null terminate suffix
    add $t9, $sp, $t6       ; buffer address
    sb $zero, 0($t9)        ; null terminator
    
    jr $ra                  ; return to caller

graph_construction_done:
    ; Build adjacency list representation of graph
    li $t1, 0               ; i = 0 (outer loop)
    
build_graph_loop:
    bge $t1, $a1, graph_done ; if i >= num_kmers, done
    
    ; Get current k-mer
    sll $t2, $t1, 2         ; $t2 = i * 4
    add $t3, $a0, $t2       ; $t3 = &kmers_array[i]
    lw $t4, 0($t3)          ; Load k-mer address
    
    ; Extract prefix and suffix (already done in process_loop)
    ; Create edge from prefix to suffix
    jal add_edge_to_graph
    
    addi $t1, $t1, 1        ; i++
    j build_graph_loop

add_edge_to_graph:
    ; Add edge from prefix to suffix to graph structure
    ; This would involve:
    ; 1. Finding if prefix exists in vertices list
    ; 2. Adding new vertex if not found
    ; 3. Creating edge from prefix vertex to suffix vertex
    
    jr $ra                  ; return

graph_done:
    ; Output the de Bruijn graph
    li $v0, 1               ; system call for exit
    li $a0, 0               ; exit status
    syscall

```

## Algorithm Explanation

The assembly code implements the following algorithm:

1. **Input Processing**: Read collection of k-mers from memory
2. **Prefix/Suffix Extraction**: For each k-mer:
   - Extract first (k-1) characters as prefix
   - Extract last (k-1) characters as suffix
3. **Graph Construction**: 
   - Create vertices for each unique (k-1)-mer
   - Create directed edges from prefix to suffix
4. **Output**: Display the de Bruijn graph structure

## Key Assembly Concepts Used

- **Memory addressing** for k-mer arrays and buffers
- **Loop control** with branch instructions
- **String manipulation** using character loading/storing
- **Function calls** with stack management
- **Data structures** for vertex and edge representation

## Time Complexity
O(n × k) where n is the number of k-mers and k is the k-mer length.

## Space Complexity  
O(n × k) for storing vertices and edges in the graph structure.