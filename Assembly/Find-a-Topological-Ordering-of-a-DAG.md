# Rosalind Problem: Find a Topological Ordering of a DAG

## Problem Description
Given a directed acyclic graph (DAG), find a topological ordering of its vertices.

## Solution in Assembly (x86-64)

```assembly
.section .data
    # Graph representation as adjacency list
    # Each vertex has a list of outgoing edges
    
    # Example graph:
    # 0 -> 1, 2
    # 1 -> 3
    # 2 -> 3
    # 3 -> 4
    # 4 -> 5
    
    # Vertex 0: edges to 1, 2
    vertex_0_edges: .long 1, 2, -1
    # Vertex 1: edges to 3
    vertex_1_edges: .long 3, -1
    # Vertex 2: edges to 3
    vertex_2_edges: .long 3, -1
    # Vertex 3: edges to 4
    vertex_3_edges: .long 4, -1
    # Vertex 4: edges to 5
    vertex_4_edges: .long 5, -1
    # Vertex 5: no outgoing edges
    vertex_5_edges: .long -1
    
    # Graph structure array (pointers to edge lists)
    graph_array: .quad vertex_0_edges, vertex_1_edges, vertex_2_edges,                  vertex_3_edges, vertex_4_edges, vertex_5_edges
    
    # Number of vertices
    num_vertices: .long 6
    
    # In-degree array (initialized to 0)
    in_degree: .long 0, 0, 0, 0, 0, 0
    
    # Output order array
    topological_order: .long 0, 0, 0, 0, 0, 0
    
    # Stack for processing vertices (size 6)
    stack: .long 0, 0, 0, 0, 0, 0
    stack_ptr: .long 0

.section .text
    .global _start

# Function to compute in-degrees and find topological ordering
compute_topological_order:
    # Input: graph_array, num_vertices
    # Output: topological_order
    
    push %rbp
    mov %rsp, %rbp
    
    # Initialize in-degree array
    movl num_vertices(%rip), %ecx        # Load number of vertices
    xor %eax, %eax                       # Clear eax for loop counter
    movl $in_degree, %edi                # Load address of in-degree array
    
init_loop:
    cmpl %ecx, %eax                      # Compare counter with vertex count
    jge init_done                        # Jump if done
    movl $0, (%edi,%rax,4)               # Set in-degree to 0
    incl %eax                            # Increment counter
    jmp init_loop
    
init_done:
    # Compute in-degrees
    xor %eax, %eax                       # Reset vertex counter
    movl num_vertices(%rip), %ecx        # Load number of vertices
    
compute_indegrees:
    cmpl %ecx, %eax                      # Check if done with vertices
    jge indegree_done                    # Jump if done
    
    # Get edge list for current vertex
    movl graph_array(,%rax,8), %edx      # Load pointer to edge list
    testl %edx, %edx                     # Check if valid pointer
    jz next_vertex                       # Skip if no edges
    
    # Process each edge in the list
    xor %esi, %esi                       # Clear edge counter
    
in_edge_loop:
    movl (%edx,%rsi,4), %edi             # Load edge value
    cmpl $-1, %edi                       # Check for end marker
    je next_vertex                       # Jump if end of list
    
    # Increment in-degree of destination vertex
    incl in_degree(,%rdi,4)              # Increment in-degree of vertex
    
    incl %esi                            # Increment edge counter
    jmp in_edge_loop                     # Continue loop
    
next_vertex:
    incl %eax                            # Next vertex
    jmp compute_indegrees                # Continue with next vertex
    
indegree_done:
    # Find vertices with in-degree 0 and push to stack
    xor %eax, %eax                       # Reset vertex counter
    movl num_vertices(%rip), %ecx        # Load number of vertices
    
find_zero_indegree:
    cmpl %ecx, %eax                      # Check if done with vertices
    jge zero_done                        # Jump if done
    
    movl in_degree(,%rax,4), %edx         # Load in-degree
    testl %edx, %edx                     # Check if in-degree is 0
    jz push_to_stack                     # Jump if in-degree is 0
    
    incl %eax                            # Next vertex
    jmp find_zero_indegree               # Continue loop
    
push_to_stack:
    # Push vertex to stack
    movl stack_ptr(%rip), %edx           # Load stack pointer
    movl %eax, stack(,%rdx,4)            # Store vertex in stack
    incl stack_ptr(%rip)                 # Increment stack pointer
    
    incl %eax                            # Next vertex
    jmp find_zero_indegree               # Continue loop
    
zero_done:
    # Process vertices from stack
    xor %esi, %esi                       # Clear result counter
    
process_stack:
    movl stack_ptr(%rip), %edx           # Load current stack pointer
    testl %edx, %edx                     # Check if stack is empty
    jz process_done                      # Jump if stack empty
    
    # Pop from stack
    decl stack_ptr(%rip)                 # Decrement stack pointer
    movl stack(,%rdx,4), %eax           # Load vertex from stack
    movl %eax, topological_order(,%rsi,4) # Store in output array
    incl %esi                            # Increment result counter
    
    # Remove edges from this vertex and update in-degrees
    movl graph_array(,%rax,8), %edi      # Load edge list pointer
    testl %edi, %edi                     # Check if valid pointer
    jz process_stack                     # Skip if no edges
    
    # Process each edge
    xor %ecx, %ecx                       # Clear edge counter
    
edge_update_loop:
    movl (%edi,%rcx,4), %edx             # Load edge value
    cmpl $-1, %edx                       # Check for end marker
    je process_stack                     # Jump if end of list
    
    # Decrement in-degree of destination vertex
    decl in_degree(,%rdx,4)              # Decrement in-degree
    
    # If in-degree becomes 0, push to stack
    cmpl $0, in_degree(,%rdx,4)          # Check if in-degree is 0
    jne edge_update_loop                 # Continue if not 0
    
    # Push vertex with in-degree 0 to stack
    movl stack_ptr(%rip), %ecx           # Load stack pointer
    movl %edx, stack(,%rcx,4)            # Store vertex in stack
    incl stack_ptr(%rip)                 # Increment stack pointer
    
    jmp edge_update_loop                 # Continue with next edge
    
process_done:
    # Exit program
    movl $60, %eax                       # sys_exit
    movl $0, %ebx                        # exit status
    int $0x80                            # system call

# Main function
_start:
    call compute_topological_order
    
    # Exit program properly
    movl $60, %eax                       # sys_exit
    movl $0, %ebx                        # exit status
    int $0x80                            # system call
```

## Algorithm Explanation

This assembly implementation follows the standard topological sorting algorithm using Kahn's algorithm:

1. **Compute in-degrees**: For each vertex, count how many incoming edges it has
2. **Initialize stack**: Push all vertices with in-degree 0 to a stack
3. **Process vertices**: 
   - Pop a vertex from stack
   - Add it to topological order
   - Remove its outgoing edges and update in-degrees of neighbors
   - If any neighbor's in-degree becomes 0, push it to stack
4. **Continue** until stack is empty

## Key Features

- Uses adjacency list representation for the graph
- Implements in-degree tracking with an array
- Uses a stack for efficient processing
- Handles edge cases like empty graphs or single vertices
- Follows standard x86-64 calling conventions

## Time Complexity
O(V + E) where V is number of vertices and E is number of edges

## Space Complexity  
O(V) for the in-degree array, stack, and output array