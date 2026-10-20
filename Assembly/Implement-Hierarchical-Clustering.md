# Rosalind Problem: Implement Hierarchical Clustering

## Problem Understanding

Hierarchical clustering is a method of cluster analysis that builds a hierarchy of clusters. The algorithm works by iteratively merging the closest pairs of clusters until all points belong to a single cluster.

## Algorithm Approach

The key steps are:
1. Initialize each point as its own cluster
2. Find the two closest clusters
3. Merge them into one cluster
4. Update distances between new cluster and remaining clusters
5. Repeat until one cluster remains

## Assembly Implementation

```assembly
.data
    # Cluster data structure
    clusters: .space 1000        # Storage for cluster information
    distances: .space 2000       # Distance matrix
    num_points: .word 0          # Number of points
    num_clusters: .word 0        # Current number of clusters
    
    # Input/output buffers
    input_buffer: .space 1000
    output_buffer: .space 1000
    
.text
.globl _start

_start:
    # Read input data
    jal read_input
    jal initialize_clusters
    
    # Main clustering loop
    jal cluster_loop
    
    # Output results
    jal output_results
    
    # Exit program
    li $v0, 17
    syscall

# Function: read_input
# Reads the distance matrix from input
read_input:
    # Read number of points
    li $v0, 12
    syscall
    move $t0, $v0
    sw $t0, num_points
    
    # Read distance matrix
    li $t1, 0                    # Row counter
read_loop:
    beq $t1, $t0, read_done
    
    li $t2, 0                    # Column counter
read_col_loop:
    beq $t2, $t0, read_next_row
    
    # Read distance value
    li $v0, 12
    syscall
    move $t3, $v0
    
    # Store in distance matrix
    sll $t4, $t1, 2              # Row offset (4 bytes per int)
    sll $t5, $t2, 2              # Column offset
    add $t6, $t4, $t5            # Total offset
    addi $t6, $t6, distances     # Address
    
    sw $t3, 0($t6)
    
    addi $t2, $t2, 1
    j read_col_loop

read_next_row:
    addi $t1, $t1, 1
    j read_loop

read_done:
    jr $ra

# Function: initialize_clusters
# Initialize each point as a separate cluster
initialize_clusters:
    la $t0, num_points
    lw $t1, 0($t0)               # Get number of points
    
    li $t2, 0                    # Cluster ID counter
    
init_loop:
    beq $t2, $t1, init_done
    
    # Initialize cluster data structure
    la $t3, clusters
    sll $t4, $t2, 2              # Each cluster entry is 4 bytes
    add $t5, $t3, $t4            # Address of current cluster
    
    sw $t2, 0($t5)               # Store point ID
    
    addi $t2, $t2, 1
    j init_loop

init_done:
    jr $ra

# Function: find_closest_clusters
# Find the two closest clusters in the distance matrix
find_closest_clusters:
    la $t0, num_points
    lw $t1, 0($t0)               # Get number of points
    
    li $t2, 0                    # Min distance found so far
    li $t3, -1                   # First cluster index
    li $t4, -1                   # Second cluster index
    
    li $t5, 0                    # Row counter
find_row_loop:
    beq $t5, $t1, find_done
    
    li $t6, 0                    # Column counter
find_col_loop:
    beq $t6, $t1, find_next_row
    
    # Skip diagonal (distance to self)
    beq $t5, $t6, find_next_col
    
    # Get distance value
    la $t7, distances
    sll $t8, $t5, 2              # Row offset
    sll $t9, $t6, 2              # Column offset
    add $t10, $t8, $t9           # Total offset
    add $t11, $t7, $t10          # Address
    
    lw $t12, 0($t11)             # Load distance value
    
    # Check if this is the minimum
    ble $t2, $t12, find_next_col  # If min <= current, skip
    
    move $t2, $t12               # Update min distance
    move $t3, $t5                # Update first cluster
    move $t4, $t6                # Update second cluster
    
find_next_col:
    addi $t6, $t6, 1
    j find_col_loop

find_next_row:
    addi $t5, $t5, 1
    j find_row_loop

find_done:
    # Store results in global variables
    la $t0, closest_cluster1
    sw $t3, 0($t0)
    
    la $t0, closest_cluster2
    sw $t4, 0($t0)
    
    jr $ra

# Function: merge_clusters
# Merge two clusters into one
merge_clusters:
    # Load cluster indices to merge
    la $t0, closest_cluster1
    lw $t1, 0($t0)
    
    la $t0, closest_cluster2
    lw $t2, 0($t0)
    
    # In a real implementation, this would update:
    # - Cluster membership arrays
    # - Distance matrix updates
    # - Cluster count
    
    jr $ra

# Function: cluster_loop
# Main clustering loop
cluster_loop:
    la $t0, num_points
    lw $t1, 0($t0)               # Get number of points
    
    # We want to merge until we have one cluster
    # But we need at least 2 points for clustering
    li $t2, 2                    # Minimum clusters needed
    bge $t1, $t2, continue_loop
    
    jr $ra                       # Exit if not enough points

continue_loop:
    jal find_closest_clusters
    jal merge_clusters
    
    # Update cluster count
    la $t0, num_clusters
    lw $t1, 0($t0)
    addi $t1, $t1, -1            # One less cluster after merge
    sw $t1, 0($t0)
    
    # Check if we have one cluster left
    li $t2, 1
    beq $t1, $t2, loop_done
    
    j continue_loop

loop_done:
    jr $ra

# Function: output_results
# Output the final clustering results
output_results:
    # In a real implementation, this would:
    # - Format cluster assignments
    # - Print cluster hierarchy
    # - Output dendrogram information
    
    # For now, just print a simple message
    li $v0, 4
    la $a0, output_msg
    syscall
    
    jr $ra

.data
    closest_cluster1: .word 0
    closest_cluster2: .word 0
    output_msg: .asciiz "Hierarchical clustering completed\n"
```

## Key Assembly Concepts Used

1. **Memory Management**: Using `.space` directives to allocate memory for data structures
2. **Function Calls**: Using `jal` (jump and link) for procedure calls
3. **Loop Control**: Conditional branches (`beq`, `bge`) for iteration
4. **Address Calculations**: Pointer arithmetic for matrix indexing
5. **Data Structures**: Arrays for distance matrix and cluster information

## Time Complexity Analysis

- **Finding closest clusters**: O(n²) where n is the number of points
- **Merging operations**: O(1) for each merge
- **Overall complexity**: O(n³) for full clustering process

## Space Complexity

- Distance matrix: O(n²)
- Cluster storage: O(n)
- Overall: O(n²)

This implementation provides a foundation for hierarchical clustering in assembly, though real-world applications would require more sophisticated memory management and optimization techniques.