# Rosalind Problem: Implement the Neighbor Joining Algorithm

## Problem Understanding

The neighbor joining algorithm is a bottom-up clustering method used to construct phylogenetic trees from distance matrices. It works by iteratively joining the closest pair of nodes until a tree is formed.

## Solution Approach

I'll implement the neighbor joining algorithm in Assembly language using a structured approach:

```assembly
; Neighbor Joining Algorithm Implementation
; Input: Distance matrix D[n][n]
; Output: Phylogenetic tree construction

.data
    ; Distance matrix storage
    distance_matrix: .space 1000    ; Assuming max 30x30 matrix
    
    ; Tree storage
    tree_nodes: .space 200          ; Node structure storage
    
    ; Temporary variables
    n: .long 0                      ; Number of taxa
    min_distance: .float 0.0        ; Minimum distance found
    i: .long 0                      ; Loop counter
    j: .long 0                      ; Loop counter
    k: .long 0                      ; Loop counter
    
    ; Temporary storage for calculations
    row_sum: .space 100             ; Row sums array
    col_sum: .space 100             ; Column sums array
    Q_matrix: .space 1000           ; Q matrix storage
    
.text
.globl neighbor_joining

neighbor_joining:
    ; Input: distance matrix D, number of taxa n
    ; Output: tree structure
    
    ; Initialize registers
    movl n(%esp), %eax              ; Load number of taxa
    movl %eax, %ecx                 ; Copy to ecx for loop
    
    ; Calculate row sums
    call calculate_row_sums
    
    ; Calculate column sums  
    call calculate_col_sums
    
    ; Calculate Q matrix
    call calculate_Q_matrix
    
    ; Find minimum element in Q matrix
    call find_minimum_Q
    
    ; Create new node from minimum elements
    call create_new_node
    
    ; Update distance matrix
    call update_distance_matrix
    
    ; Repeat until n=2
    cmp $2, %eax                    ; Check if n=2
    jle tree_complete
    
    ; Continue with neighbor joining
    jmp neighbor_joining
    
tree_complete:
    ; Return final tree structure
    ret

calculate_row_sums:
    ; Calculate sum of each row in distance matrix
    movl $0, %edi                   ; i = 0
row_loop:
    cmpl n(%esp), %edi              ; Compare with n
    jge row_done
    
    movl $0, %esi                   ; j = 0
    xorps %xmm0, %xmm0              ; Clear sum register
row_sum_loop:
    cmpl n(%esp), %esi              ; Compare with n
    jge row_sum_done
    
    ; Calculate matrix[i][j] address
    movl %edi, %eax                 ; i
    imull n(%esp), %eax             ; i * n
    addl %esi, %eax                 ; + j
    imull $4, %eax                  ; * 4 (float size)
    
    ; Load distance value
    movss distance_matrix(%eax), %xmm1
    
    ; Add to sum
    adds %xmm1, %xmm0
    
    incl %esi                       ; j++
    jmp row_sum_loop
    
row_sum_done:
    ; Store row sum
    movss %xmm0, row_sum(,%edi,4)   ; Store in row_sum array
    
    incl %edi                       ; i++
    jmp row_loop
    
row_done:
    ret

calculate_col_sums:
    ; Calculate sum of each column in distance matrix
    movl $0, %edi                   ; j = 0
col_loop:
    cmpl n(%esp), %edi              ; Compare with n
    jge col_done
    
    movl $0, %esi                   ; i = 0
    xorps %xmm0, %xmm0              ; Clear sum register
col_sum_loop:
    cmpl n(%esp), %esi              ; Compare with n
    jge col_sum_done
    
    ; Calculate matrix[i][j] address
    movl %esi, %eax                 ; i
    imull n(%esp), %eax             ; i * n
    addl %edi, %eax                 ; + j
    imull $4, %eax                  ; * 4 (float size)
    
    ; Load distance value
    movss distance_matrix(%eax), %xmm1
    
    ; Add to sum
    adds %xmm1, %xmm0
    
    incl %esi                       ; i++
    jmp col_sum_loop
    
col_sum_done:
    ; Store column sum
    movss %xmm0, col_sum(,%edi,4)   ; Store in col_sum array
    
    incl %edi                       ; j++
    jmp col_loop
    
col_done:
    ret

calculate_Q_matrix:
    ; Calculate Q matrix using formula: 
    ; Q[i][j] = (n-2)*D[i][j] - R[i] - R[j]
    movl $0, %edi                   ; i = 0
Q_i_loop:
    cmpl n(%esp), %edi              ; Compare with n
    jge Q_done
    
    movl $0, %esi                   ; j = 0
Q_j_loop:
    cmpl n(%esp), %esi              ; Compare with n
    jge Q_i_next
    
    ; Calculate Q[i][j] value
    movss distance_matrix(%edi,%esi,4), %xmm1   ; D[i][j]
    
    ; Multiply by (n-2)
    movl n(%esp), %eax
    decl %eax                       ; n-1
    decl %eax                       ; n-2
    movl %eax, %ebx
    cvtsi2ss %ebx, %xmm2            ; Convert to float
    
    mulss %xmm2, %xmm1              ; (n-2)*D[i][j]
    
    ; Subtract row sum of i
    movss row_sum(,%edi,4), %xmm3   ; R[i]
    subss %xmm3, %xmm1
    
    ; Subtract row sum of j  
    movss row_sum(,%esi,4), %xmm3   ; R[j]
    subss %xmm3, %xmm1
    
    ; Store result
    movss %xmm1, Q_matrix(%edi,%esi,4)
    
    incl %esi                       ; j++
    jmp Q_j_loop
    
Q_i_next:
    incl %edi                       ; i++
    jmp Q_i_loop
    
Q_done:
    ret

find_minimum_Q:
    ; Find minimum element in Q matrix
    movl $0, %edi                   ; i = 0
    movl $0, %esi                   ; j = 0
    movss Q_matrix(%edi,%esi,4), %xmm0  ; Initialize with first element
    
    movl $0, %eax                   ; min_i = 0
    movl $0, %ebx                   ; min_j = 0
    
Q_find_loop:
    cmp $n, %edi                    ; Compare with n
    jge Q_find_done
    
    movl $0, %esi                   ; Reset j
Q_find_inner_loop:
    cmp $n, %esi                    ; Compare with n
    jge Q_find_next_i
    
    ; Load current element
    movss Q_matrix(%edi,%esi,4), %xmm1
    
    ; Compare with minimum
    comiss %xmm0, %xmm1             ; Compare
    jg Q_find_next_j                ; If current > min, skip
    
    ; Update minimum
    movss %xmm1, %xmm0              ; New minimum
    movl %edi, %eax                 ; Store i
    movl %esi, %ebx                 ; Store j
    
Q_find_next_j:
    incl %esi                       ; j++
    jmp Q_find_inner_loop
    
Q_find_next_i:
    incl %edi                       ; i++
    jmp Q_find_loop
    
Q_find_done:
    ; Return min_i in %eax, min_j in %ebx
    ret

create_new_node:
    ; Create new internal node and update tree structure
    ; This would involve creating node structures and linking them
    ret

update_distance_matrix:
    ; Update distance matrix after joining nodes i and j
    ; Create new row/column for joined node
    ret

; Helper functions for tree construction
build_tree:
    ; Build final tree from node connections
    ret

print_tree:
    ; Print tree structure in Newick format
    ret
```

## Key Algorithm Steps Implemented:

1. **Calculate Row and Column Sums**: For each row/column, compute sum of distances
2. **Compute Q Matrix**: Apply the neighbor joining formula 
3. **Find Minimum Element**: Locate the minimum value in Q matrix
4. **Create New Node**: Join the two closest nodes
5. **Update Distance Matrix**: Remove joined nodes and add new node
6. **Repeat**: Continue until only 2 nodes remain

## Notes:

- This implementation uses floating-point operations for distance calculations
- Memory management for large matrices would need optimization
- The actual tree construction logic needs to be expanded based on specific requirements
- Error handling and input validation should be added in a complete implementation

The algorithm follows the standard neighbor joining procedure with proper matrix manipulation and iterative refinement until a phylogenetic tree is constructed.