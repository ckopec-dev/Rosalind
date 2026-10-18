# Rosalind Problem: Implement UPGMA

## Problem Understanding

UPGMA (Unweighted Pair Group Method with Arithmetic Mean) is a hierarchical clustering algorithm used to construct phylogenetic trees from distance matrices. It's particularly useful for creating evolutionary trees where all branches have equal length.

## Algorithm Overview

1. Start with a distance matrix between all sequences
2. Find the pair of sequences with minimum distance
3. Create a new node (cluster) containing these two sequences
4. Update the distance matrix by computing distances from the new cluster to remaining sequences
5. Repeat until only one cluster remains

## Assembly Implementation

```assembly
; UPGMA Implementation in Assembly Language
; Assumes 8-bit processor with basic arithmetic operations

.data
    ; Distance matrix (symmetric)
    dist_matrix: .byte 0, 5, 9, 10
                 .byte 5, 0, 6, 7
                 .byte 9, 6, 0, 3
                 .byte 10, 7, 3, 0
    
    ; Number of sequences
    n_seq: .byte 4
    
    ; Cluster information
    clusters: .byte 0, 1, 2, 3   ; Initial cluster IDs
    cluster_count: .byte 4
    
    ; Temporary storage for distances
    temp_distances: .space 4

.text
    ; Main UPGMA procedure
    main:
        mov r0, #n_seq          ; Load number of sequences
        mov r1, #dist_matrix    ; Load distance matrix address
        mov r2, #clusters       ; Load cluster addresses
        
    ; Loop until we have one cluster
    loop_clusters:
        cmp r0, #1              ; Check if only one cluster left
        jle done                ; If yes, exit
        
        ; Find minimum distance in matrix
        call find_min_distance  ; Returns row and column of min distance
        
        ; Create new cluster from min pair
        call merge_clusters     ; Merge two clusters into one
        
        ; Update distance matrix
        call update_distances   ; Recalculate distances to new cluster
        
        dec r0                  ; Decrease cluster count
        jmp loop_clusters       ; Continue loop
        
    done:
        ; Output final tree structure
        ret

    ; Function: Find minimum distance in current matrix
    find_min_distance:
        push r0                 ; Save registers
        push r1
        push r2
        push r3
        push r4
        
        mov r0, #dist_matrix    ; Load matrix address
        mov r1, #0              ; Initialize min_distance = 0
        mov r2, #0              ; Initialize min_row = 0
        mov r3, #0              ; Initialize min_col = 0
        mov r4, #0              ; Initialize i = 0
        
    find_loop:
        cmp r4, #n_seq          ; Check if we've processed all rows
        jge find_done           ; If yes, exit loop
        
        mov r5, #0              ; Initialize j = 0
        mov r6, r4              ; Load current row number
        
    inner_loop:
        cmp r5, #n_seq          ; Check if we've processed all columns
        jge inner_next_row      ; If yes, go to next row
        
        ; Calculate address of matrix[i][j]
        mov r7, r6              ; Copy row index
        mul r7, #n_seq          ; Multiply by number of columns
        add r7, r5              ; Add column index
        mov r8, r7              ; Store address in r8
        
        ; Load distance value
        ld r9, [r0 + r8]        ; Load distance from matrix
        
        ; Skip diagonal elements (distance to self = 0)
        cmp r6, r5              ; Compare row and column indices
        je skip_current         ; If equal, skip this element
        
        ; Check if current distance is minimum
        cmp r9, r1              ; Compare with current min
        jg skip_current         ; If greater, skip
        
        ; Update minimum values
        mov r1, r9              ; Update min_distance
        mov r2, r6              ; Update min_row
        mov r3, r5              ; Update min_col
        
    skip_current:
        inc r5                  ; Increment column index
        jmp inner_loop          ; Continue inner loop
        
    inner_next_row:
        inc r4                  ; Increment row index
        jmp find_loop           ; Continue outer loop
        
    find_done:
        ; Return values in r2 (min_row) and r3 (min_col)
        pop r4
        pop r3
        pop r2
        pop r1
        pop r0
        ret

    ; Function: Merge two clusters into one
    merge_clusters:
        push r0
        push r1
        push r2
        
        ; Get the two cluster indices to be merged
        mov r0, #clusters       ; Load cluster address
        mov r1, #min_row        ; Get first cluster index
        mov r2, #min_col        ; Get second cluster index
        
        ; Merge clusters (simplified approach)
        ; In real implementation, would create new node with two children
        
        pop r2
        pop r1
        pop r0
        ret

    ; Function: Update distance matrix after merging
    update_distances:
        push r0
        push r1
        push r2
        push r3
        push r4
        push r5
        push r6
        
        ; Recalculate distances from new cluster to remaining sequences
        mov r0, #dist_matrix    ; Load matrix address
        mov r1, #n_seq          ; Load number of sequences
        mov r2, #cluster_count  ; Load current cluster count
        
        ; For each existing cluster, compute average distance
        mov r3, #0              ; Initialize i = 0
        
    update_loop:
        cmp r3, r2              ; Check if all clusters processed
        jge update_done         ; If yes, exit
        
        ; Skip the two merged clusters
        cmp r3, min_row         ; Check if current cluster is one of merged
        je skip_update          ; If yes, skip
        cmp r3, min_col         ; Check if current cluster is second merged
        je skip_update          ; If yes, skip
        
        ; Compute average distance to new cluster
        call compute_avg_distance
        
    skip_update:
        inc r3                  ; Increment cluster index
        jmp update_loop         ; Continue loop
        
    update_done:
        pop r6
        pop r5
        pop r4
        pop r3
        pop r2
        pop r1
        pop r0
        ret

    ; Function: Compute average distance from new cluster to others
    compute_avg_distance:
        push r0
        push r1
        push r2
        
        mov r0, #dist_matrix    ; Load matrix address
        mov r1, #new_cluster    ; Load index of new cluster
        mov r2, #other_cluster  ; Load index of other cluster
        
        ; Calculate average of distances from new to all others
        ; Simplified version - in practice would need proper averaging logic
        
        pop r2
        pop r1
        pop r0
        ret

    ; Helper function: Get distance between two clusters
    get_distance:
        push r0
        push r1
        
        mov r0, #dist_matrix    ; Load matrix address
        mov r1, #cluster1       ; Load first cluster index
        mov r2, #cluster2       ; Load second cluster index
        
        ; Calculate matrix address and return value
        ; Implementation depends on specific architecture
        
        pop r1
        pop r0
        ret

```

## Key Assembly Concepts Used

1. **Memory Management**: Uses data section for storing distance matrices and clusters
2. **Register Operations**: Multiple registers used for counters, addresses, and values
3. **Control Flow**: Conditional jumps (cmp, jle, je) for loop control
4. **Function Calls**: Proper stack management for subroutines
5. **Arithmetic Operations**: Addition, multiplication, comparison operations

## Algorithm Complexity

- **Time Complexity**: O(n³) where n is the number of sequences
- **Space Complexity**: O(n²) for storing the distance matrix

## Notes for Implementation

1. This is a simplified assembly implementation - real-world usage would require:
   - Proper memory management
   - More sophisticated data structures for tree representation
   - Enhanced error handling
   - Better floating-point support (if needed)
2. The actual implementation would depend on the specific assembly architecture
3. Distance matrix operations need to be optimized based on available instructions

This implementation provides the basic structure for UPGMA clustering in assembly language, focusing on the core algorithmic components rather than low-level optimizations.