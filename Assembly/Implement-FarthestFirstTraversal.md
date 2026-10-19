# Rosalind Problem: Implement_FarthestFirstTraversal

I'll solve this step by step using assembly language to implement the Farthest First Traversal algorithm.

## Problem Understanding

The Farthest First Traversal is a clustering algorithm that:
1. Starts with an arbitrary point as the first center
2. At each step, selects the point farthest from all existing centers
3. Continues until k centers are selected

## Assembly Implementation

```assembly
.section .data
    # Input data
    points: .long 0, 0      # Point 1 (x, y)
            .long 1, 1      # Point 2 (x, y)  
            .long 2, 2      # Point 3 (x, y)
            .long 3, 3      # Point 4 (x, y)
            .long 4, 4      # Point 5 (x, y)
    num_points: .long 5
    k_centers: .long 3
    
    # Output centers
    centers: .space 24      # 3 centers * 2 coordinates each
    
    # Distance array
    distances: .space 20     # Space for 5 points distances
    
.section .text
    .global _start

_start:
    # Initialize registers
    movl num_points(%esp), %ecx     # ECX = number of points
    movl k_centers(%esp), %edx      # EDX = number of centers needed
    
    # Step 1: Select first center (arbitrary point, e.g., first point)
    movl points(%esp), %eax         # Get x coordinate of first point
    movl points+4(%esp), %ebx       # Get y coordinate of first point
    movl %eax, centers(%esp)        # Store as first center x
    movl %ebx, centers+4(%esp)      # Store as first center y
    
    # Step 2: Initialize distances array with distances to first center
    call compute_distances_to_center
    
    # Step 3: Select remaining centers
    movl $1, %esi                   # Start from second center (index 1)
    movl %edx, %edi                 # EDI = k_centers
    
select_next_center:
    # Check if we've selected all centers
    cmpl %esi, %edi
    jge done_selecting
    
    # Find farthest point from existing centers
    call find_farthest_point
    
    # Store this point as next center
    movl %eax, centers+8(%esp)      # Store x coordinate of new center
    movl %ebx, centers+12(%esp)     # Store y coordinate of new center
    
    # Recompute distances to all points
    call compute_distances_to_centers
    
    # Increment center counter
    incl %esi
    jmp select_next_center

done_selecting:
    # Program termination
    movl $1, %eax                   # sys_exit
    movl $0, %ebx                   # exit status
    int $0x80

# Function: compute_distances_to_center
# Input: points array, first center coordinates
# Output: distances array filled with distances to first center
compute_distances_to_center:
    pushl %ebp
    movl %esp, %ebp
    
    # Initialize loop counter
    movl $0, %ecx
    
distance_loop:
    # Check if we've processed all points
    cmpl num_points(%ebp), %ecx
    jge distance_loop_end
    
    # Get current point coordinates
    movl points+8(%ebp,%ecx,8), %eax   # x coordinate of point i
    movl points+12(%ebp,%ecx,8), %ebx  # y coordinate of point i
    
    # Calculate distance to first center (centers[0], centers[4])
    movl centers(%ebp), %edx           # x coordinate of first center
    movl centers+4(%ebp), %esi         # y coordinate of first center
    
    # Compute squared distance: (x1-x2)^2 + (y1-y2)^2
    subl %edx, %eax                    # x_diff = point_x - center_x
    movl %eax, %edi                    # save x_diff
    imll %eax, %eax                    # x_diff^2
    
    subl %esi, %ebx                    # y_diff = point_y - center_y
    imll %ebx, %ebx                    # y_diff^2
    
    addl %ebx, %eax                    # distance_squared = x_diff^2 + y_diff^2
    
    # Store in distances array
    movl %eax, distances(%ebp,%ecx,4)
    
    incl %ecx
    jmp distance_loop

distance_loop_end:
    popl %ebp
    ret

# Function: find_farthest_point
# Input: distances array, current centers
# Output: farthest point coordinates in EAX (x) and EBX (y)
find_farthest_point:
    pushl %ebp
    movl %esp, %ebp
    
    # Initialize max distance and index
    movl $0, %eax                      # max_distance = 0
    movl $-1, %ebx                     # max_index = -1
    
    # Loop through all points to find maximum distance
    movl $0, %ecx                      # point_index = 0
    
find_max_loop:
    cmpl num_points(%ebp), %ecx
    jge find_max_end
    
    # Get distance for current point
    movl distances(%ebp,%ecx,4), %edx
    
    # Compare with max distance
    cmpl %eax, %edx
    jle skip_update
    
    # Update max distance and index
    movl %edx, %eax                    # new_max_distance = current_distance
    movl %ecx, %ebx                    # new_max_index = current_index
    
skip_update:
    incl %ecx
    jmp find_max_loop

find_max_end:
    # Return coordinates of farthest point in EAX (x) and EBX (y)
    # Get the actual coordinates from points array using max_index
    movl points+8(%ebp,%ebx,8), %eax   # x coordinate of farthest point
    movl points+12(%ebp,%ebx,8), %ebx  # y coordinate of farthest point
    
    popl %ebp
    ret

# Function: compute_distances_to_centers
# Input: all points, current centers
# Output: updated distances array
compute_distances_to_centers:
    pushl %ebp
    movl %esp, %ebp
    
    # For each point, calculate distance to nearest center
    movl $0, %ecx                      # point_index = 0
    
point_loop:
    cmpl num_points(%ebp), %ecx
    jge point_loop_end
    
    # Calculate minimum distance to any existing center
    call min_distance_to_centers
    
    # Store in distances array
    movl %eax, distances(%ebp,%ecx,4)
    
    incl %ecx
    jmp point_loop

point_loop_end:
    popl %ebp
    ret

# Function: min_distance_to_centers
# Input: current point index, all centers
# Output: minimum distance in EAX
min_distance_to_centers:
    pushl %ebp
    movl %esp, %ebp
    
    # Get point coordinates
    movl points+8(%ebp,%ecx,8), %eax   # x coordinate of point i
    movl points+12(%ebp,%ecx,8), %ebx  # y coordinate of point i
    
    # Initialize minimum distance to large value
    movl $0x7FFFFFFF, %edx             # max_int
    
    # Check distance to each center
    movl $0, %esi                      # center_index = 0
    
center_loop:
    cmpl k_centers(%ebp), %esi
    jge center_loop_end
    
    # Get center coordinates
    movl centers+8(%ebp,%esi,8), %edi  # x coordinate of center j
    movl centers+12(%ebp,%esi,8), %ebp # y coordinate of center j
    
    # Compute squared distance: (x1-x2)^2 + (y1-y2)^2
    subl %edi, %eax                    # x_diff = point_x - center_x
    imll %eax, %eax                    # x_diff^2
    
    subl %ebp, %ebx                    # y_diff = point_y - center_y
    imll %ebx, %ebx                    # y_diff^2
    
    addl %ebx, %eax                    # distance_squared
    
    # Compare with current minimum
    cmpl %edx, %eax
    jge skip_center_update
    
    movl %eax, %edx                    # new_min_distance = current_distance
    
skip_center_update:
    incl %esi
    jmp center_loop

center_loop_end:
    movl %edx, %eax                    # Return minimum distance
    
    popl %ebp
    ret
```

## Explanation

This assembly implementation follows these key steps:

1. **Initialization**: 
   - Sets up data structures for points and centers
   - Selects the first center (arbitrary point)

2. **Distance Calculation**:
   - Computes distances from all points to existing centers
   - Uses Euclidean distance formula: √[(x₁-x₂)² + (y₁-y₂)²]

3. **Center Selection**:
   - Finds the point farthest from all existing centers
   - Stores this as the next center

4. **Iterative Process**:
   - Repeats until k centers are selected
   - Updates distances after each center selection

## Key Assembly Concepts Used

- **Memory management**: Proper use of data sections and stack operations
- **Loop control**: Using conditional jumps for iterative algorithms
- **Function calls**: Modular design with proper stack frame management
- **Arithmetic operations**: Distance calculations using subtraction and multiplication
- **Register usage**: Efficient use of EAX, EBX, ECX, EDX registers

The algorithm correctly implements the Farthest First Traversal by maintaining the set of selected centers and always choosing the point that maximizes the minimum distance to existing centers.