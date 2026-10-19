# Rosalind Problem: Compute the Squared Error Distortion

## Problem Understanding

The squared error distortion is a measure of how well a set of data points fits a set of cluster centers (centroids). For each data point, we calculate the squared Euclidean distance to its nearest centroid and sum these distances.

## Solution in Assembly Language

```assembly
.section .data
    # Sample data points (2D coordinates)
    points: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
    num_points: .long 4
    
    # Cluster centers (centroids)
    centroids: .float 2.0, 3.0, 6.0, 7.0
    num_centroids: .long 2
    
    # Result variable
    distortion: .float 0.0

.section .text
    .global _start

_start:
    # Initialize registers
    la $t0, points          # Load address of points array
    li $t1, 0               # Point counter (i = 0)
    li $t2, 0               # Total distortion sum
    
    # Get number of points
    la $t3, num_points
    lw $t4, 0($t3)          # Load number of points
    
    # Main loop over all points
point_loop:
    # Check if we've processed all points
    bge $t1, $t4, end_calculation
    
    # Calculate squared distance to nearest centroid for current point
    jal calculate_distortion_for_point
    
    # Add distortion to total
    add $t2, $t2, $v0       # Accumulate distortion
    
    # Move to next point
    addi $t1, $t1, 1        # i++
    j point_loop

end_calculation:
    # Store final result
    la $t3, distortion
    sw $t2, 0($t3)
    
    # Exit program
    li $v0, 10              # System call for exit
    syscall

# Function: calculate_distortion_for_point
# Input: $t1 = current point index
# Output: $v0 = distortion value for this point
calculate_distortion_for_point:
    # Save registers
    push $ra
    push $s0
    push $s1
    
    # Calculate address of current point (2 floats per point)
    sll $t5, $t1, 2         # Multiply index by 4 bytes (float size)
    add $t5, $t5, $t0       # Add base address of points
    
    # Load current point coordinates
    lwc1 $f0, 0($t5)        # x coordinate of current point
    lwc1 $f2, 4($t5)        # y coordinate of current point
    
    # Initialize minimum distance and centroid index
    li $t6, 0               # centroid counter
    li $t7, 0               # min_distance (set to 0 initially)
    li $t8, -1              # best_centroid_index (-1 means not set)
    
    # Get number of centroids
    la $t9, num_centroids
    lw $t10, 0($t9)         # Load number of centroids
    
centroid_loop:
    # Check if we've processed all centroids
    bge $t6, $t10, centroid_end
    
    # Calculate address of current centroid (2 floats per centroid)
    sll $t11, $t6, 3        # Multiply index by 8 bytes (2 * 4 bytes)
    la $t12, centroids
    add $t11, $t11, $t12
    
    # Load centroid coordinates
    lwc1 $f4, 0($t11)       # x coordinate of centroid
    lwc1 $f6, 4($t11)       # y coordinate of centroid
    
    # Calculate squared Euclidean distance
    # Distance = (x1-x2)^2 + (y1-y2)^2
    sub.s $f8, $f0, $f4     # x_diff = point_x - centroid_x
    sub.s $f10, $f2, $f6    # y_diff = point_y - centroid_y
    
    mul.s $f12, $f8, $f8    # x_diff^2
    mul.s $f14, $f10, $f10  # y_diff^2
    
    add.s $f16, $f12, $f14  # squared_distance = x_diff^2 + y_diff^2
    
    # Check if this is the minimum distance so far
    li $t13, 1              # Flag for first centroid
    beq $t8, -1, set_min    # If no minimum yet, set current as minimum
    
    c.lt.s $f16, $f7        # Compare with current minimum
    bc1f not_better         # If current is not better, skip
    
set_min:
    mov.s $f7, $f16         # Set new minimum distance
    li $t8, $t6             # Store centroid index
    
not_better:
    addi $t6, $t6, 1        # Next centroid
    j centroid_loop

centroid_end:
    # Return minimum distance as distortion for this point
    mov.s $v0, $f7         # Move result to return register
    
    # Restore registers
    pop $s1
    pop $s0
    pop $ra
    jr $ra

# Helper function to push registers (simplified version)
push:
    addi $sp, $sp, -4
    sw $a0, 0($sp)
    jr $ra

# Helper function to pop registers (simplified version)
pop:
    lw $a0, 0($sp)
    addi $sp, $sp, 4
    jr $ra
```

## Alternative Implementation (More Compact)

```assembly
.section .data
    # Sample data points: [x1,y1,x2,y2,x3,y3,x4,y4]
    points: .float 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0
    num_points: .long 4
    
    # Centroids: [cx1,cy1,cx2,cy2]
    centroids: .float 2.0, 3.0, 6.0, 7.0
    num_centroids: .long 2

.section .text
    .global _start

_start:
    # Initialize variables
    la $s0, points          # Points array address
    la $s1, centroids       # Centroids array address
    li $s2, 0               # Point counter
    li $s3, 0               # Total distortion
    
    # Load number of points
    lw $s4, num_points      # Number of points
    
point_loop:
    # Check if done with all points
    bge $s2, $s4, print_result
    
    # Calculate distortion for current point
    jal calculate_squared_error
    
    # Add to total
    add $s3, $s3, $v0       # Accumulate distortion
    
    # Next point
    addi $s2, $s2, 1
    j point_loop

print_result:
    # Store final result in memory
    la $t0, distortion_result
    sw $s3, 0($t0)
    
    # Exit program
    li $v0, 10              # exit syscall
    syscall

# Function to calculate squared error for one point
calculate_squared_error:
    # Calculate address of current point (2 floats = 8 bytes)
    sll $t0, $s2, 3         # index * 8
    add $t0, $t0, $s0       # + base address
    
    # Load point coordinates
    lwc1 $f0, 0($t0)        # x coordinate
    lwc1 $f2, 4($t0)        # y coordinate
    
    # Initialize minimum distance tracking
    li $t1, 0               # centroid counter
    li $t2, 1               # flag for first centroid
    
    # Load number of centroids
    lw $t3, num_centroids
    
min_distance_loop:
    bge $t1, $t3, end_calculation
    
    # Calculate centroid address
    sll $t4, $t1, 3         # index * 8
    add $t4, $t4, $s1       # + base address
    
    # Load centroid coordinates
    lwc1 $f4, 0($t4)        # cx
    lwc1 $f6, 4($t4)        # cy
    
    # Calculate squared Euclidean distance
    sub.s $f8, $f0, $f4     # x_diff = point_x - centroid_x
    sub.s $f10, $f2, $f6    # y_diff = point_y - centroid_y
    
    mul.s $f12, $f8, $f8    # x_diff^2
    mul.s $f14, $f10, $f10  # y_diff^2
    
    add.s $f16, $f12, $f14  # distance^2
    
    # Initialize or update minimum
    beq $t2, 1, set_initial # First centroid case
    
    c.lt.s $f16, $f18       # Compare with current min
    bc1f skip_update        # If not better, skip
    
set_initial:
    mov.s $f18, $f16        # Set initial minimum
    li $t2, 0               # Not first anymore
    
skip_update:
    addi $t1, $t1, 1        # Next centroid
    j min_distance_loop

end_calculation:
    # Return minimum distance as distortion value
    mov.s $v0, $f18         # Return result in $v0
    
    jr $ra
```

## Key Concepts Used

1. **Memory Access**: Loading data points and centroids from memory
2. **Floating Point Operations**: Calculating Euclidean distances using floating-point arithmetic
3. **Loop Control**: Processing all data points and centroids
4. **Function Calls**: Modular approach to calculating distortion for each point
5. **Register Management**: Proper saving/restoring of registers during function calls

## Algorithm Steps

1. For each data point, find the nearest centroid
2. Calculate squared Euclidean distance between point and centroid
3. Sum all squared distances across all points
4. Return the total distortion value

This implementation handles the core computational requirements of the squared error distortion calculation in assembly language.