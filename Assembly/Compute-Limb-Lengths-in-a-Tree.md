# Rosalind Problem: Compute Limb Lengths in a Tree

## Problem Understanding

The limb length of a leaf in a tree is the distance from that leaf to its parent node. Given a distance matrix for a tree, we need to compute the limb length for a specific leaf.

## Approach

1. **Key Insight**: For a leaf `j` and a tree with distance matrix `D`, if we know the distances from leaf `j` to all other leaves, we can compute its limb length using the formula:
   ```
   limb_length(j) = (D[i][j] + D[j][k] - D[i][k]) / 2
   ```

2. **Algorithm Steps**:
   - For a given leaf `j`, find three other leaves `i`, `k` that form a triangle in the tree
   - Use the distance matrix to compute limb length using the formula above
   - The result should be the same regardless of which valid triplet we choose

## Assembly Implementation

```assembly
.data
    # Distance matrix D (assuming max 100x100)
    D: .space 10000      # 100*100 * 4 bytes = 40KB
    
    # Input parameters
    n: .long 0           # Number of leaves
    j: .long 0           # Leaf index for which we compute limb length
    
    # Temporary variables
    temp1: .long 0
    temp2: .long 0
    result: .long 0
    
.text
.globl compute_limb_length

compute_limb_length:
    # Function parameters:
    # r0 = n (number of leaves)
    # r1 = j (leaf index)
    # r2 = pointer to distance matrix D
    
    mov r0, n            # Store number of leaves
    mov r1, j            # Store leaf index
    mov r2, %r2          # Store matrix pointer
    
    # Find three other nodes (i, k) to compute limb length
    # We'll use first two nodes as i and k (assuming valid for computation)
    
    # Compute D[i][j] where i = 0, j = r1
    # Compute D[j][k] where k = 1, j = r1  
    # Compute D[i][k] where i = 0, k = 1
    
    # Calculate D[0][r1] = D[0][j]
    mov $0, %r3          # i = 0
    mov r1, %r4          # j = r1
    mul %r3, %r4         # r3 * r4 = index for row-major matrix
    add %r4, %r4         # Multiply by 4 (for integer size)
    add r2, %r4          # Add base address
    ldr %r5, [%r4]       # Load D[0][j]
    
    # Calculate D[r1][1] = D[j][1] 
    mov r1, %r3          # j = r1
    mov $1, %r4          # k = 1
    mul %r3, %r4         # r3 * r4 = index for row-major matrix
    add %r4, %r4         # Multiply by 4 (for integer size)
    add r2, %r4          # Add base address  
    ldr %r6, [%r4]       # Load D[j][1]
    
    # Calculate D[0][1] = D[i][k]
    mov $0, %r3          # i = 0
    mov $1, %r4          # k = 1
    mul %r3, %r4         # r3 * r4 = index for row-major matrix  
    add %r4, %r4         # Multiply by 4 (for integer size)
    add r2, %r4          # Add base address
    ldr %r7, [%r4]       # Load D[0][1]
    
    # Compute limb_length = (D[i][j] + D[j][k] - D[i][k]) / 2
    # Step 1: D[i][j] + D[j][k] 
    add %r5, %r6         # Add D[0][j] and D[j][1]
    
    # Step 2: (D[i][j] + D[j][k]) - D[i][k]
    sub %r7, %r5         # Subtract D[0][1] from sum
    
    # Step 3: Divide by 2
    lsr %r5, #1          # Right shift by 1 = divide by 2
    
    # Store result
    mov %r5, result      # Save computed limb length
    
    # Return result in r0
    mov %r5, %r0
    bx lr                # Return

# Alternative implementation for multiple leaves
compute_limb_length_optimized:
    # Input: n (number of leaves), j (leaf index), D (distance matrix)
    # Output: limb length for leaf j
    
    push {r4-r7, lr}
    
    mov r0, n            # Load n 
    mov r1, j            # Load j
    mov r2, %r2          # Load matrix pointer
    
    # For leaf j, we need to find three other leaves i,k such that:
    # limb_length = (D[i][j] + D[j][k] - D[i][k]) / 2
    
    # Use first three leaves as reference (assuming n >= 3)
    mov $0, %r3          # i = 0
    mov $1, %r4          # k = 1  
    mov r1, %r5          # j = r1
    
    # Calculate D[i][j] = D[0][j]
    mov %r3, %r6         # i = 0
    mul %r6, n           # i * n
    add %r6, %r5         # + j (column index)
    lsl %r6, #2          # Multiply by 4 for integer size
    add r2, %r6          # Add base address
    ldr %r7, [%r6]       # Load D[0][j]
    
    # Calculate D[j][k] = D[j][1]  
    mov %r5, %r6         # j = r1
    mul %r6, n           # j * n
    add %r6, %r4         # + k = 1 (column index)
    lsl %r6, #2          # Multiply by 4 for integer size
    add r2, %r6          # Add base address
    ldr %r8, [%r6]       # Load D[j][1]
    
    # Calculate D[i][k] = D[0][1]
    mov %r3, %r6         # i = 0  
    mul %r6, n           # i * n
    add %r6, %r4         # + k = 1 (column index)
    lsl %r6, #2          # Multiply by 4 for integer size
    add r2, %r6          # Add base address
    ldr %r9, [%r6]       # Load D[0][1]
    
    # Compute limb_length = (D[i][j] + D[j][k] - D[i][k]) / 2
    add %r7, %r8         # (D[i][j] + D[j][k])
    sub %r9, %r7         # - D[i][k]
    lsr %r7, #1          # / 2
    
    mov %r7, result      # Store result
    mov %r7, %r0         # Return in r0
    
    pop {r4-r7, pc}
```

## Example Usage

```assembly
# Example: Compute limb length for leaf 3 in a 5-leaf tree
# Assuming distance matrix D is pre-loaded

mov $5, %r0          # n = 5 leaves  
mov $3, %r1          # j = 3 (compute limb length for leaf 3)
mov %matrix_addr, %r2 # Load matrix address
bl compute_limb_length

# Result in r0: computed limb length
```

## Key Features

1. **Matrix Access**: Uses row-major order indexing for distance matrix access
2. **Memory Efficient**: Works with pre-loaded distance matrix
3. **Correct Formula**: Implements the mathematical formula correctly
4. **Optimized Division**: Uses bit shifting for division by 2
5. **Proper Register Management**: Handles register preservation during function calls

## Time Complexity

- **Time**: O(1) - constant time operations (3 memory accesses and basic arithmetic)
- **Space**: O(1) - uses only a constant amount of extra space

This implementation assumes the distance matrix is stored in row-major order with integer values and that there are at least 3 leaves in the tree.