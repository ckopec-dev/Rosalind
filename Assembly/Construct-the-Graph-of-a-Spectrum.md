# Construct the Graph of a Spectrum - Rosalind Problem Solution

## Problem Understanding

The task is to construct a directed graph from a spectrum of masses, where:
- Each node represents a mass value
- Edges represent transitions between masses (differences)
- We need to build a graph where each edge weight equals the difference between two masses

## Solution Approach

We'll implement this in Assembly language using a simple approach:
1. Parse input spectrum values
2. Calculate all possible differences between masses
3. Build adjacency list representation of the graph

## Assembly Implementation

```assembly
.section .data
    # Input spectrum array (example values)
    spectrum: .long 57, 71, 101, 113, 129, 147, 156, 170, 186, 194, 222, 231, 248, 263, 271, 284, 295, 305, 307, 320, 332, 341, 347, 357, 364, 369, 376, 383, 388, 391, 396, 402, 407, 413, 415, 421, 428, 430, 435, 437, 443, 447, 451, 453, 456, 459, 460, 463, 466, 468, 470, 472, 473, 475, 477, 478, 480, 482, 484, 486, 488, 490, 491, 493, 495, 497, 499, 500, 502, 504, 506, 508, 510, 512, 514, 516, 518, 520, 522, 524, 526, 528, 530, 532, 534, 536, 538, 540, 542, 544, 546, 548, 550, 552, 554, 556, 558, 560, 562, 564, 566, 568, 570, 572, 574, 576, 578, 580, 582, 584, 586, 588, 590, 592, 594, 596, 598, 600, 602, 604, 606, 608, 610, 612, 614, 616, 618, 620, 622, 624, 626, 628, 630, 632, 634, 636, 638, 640, 642, 644, 646, 648, 650, 652, 654, 656, 658, 660, 662, 664, 666, 668, 670, 672, 674, 676, 678, 680, 682, 684, 686, 688, 690, 692, 694, 696, 698, 700, 702, 704, 706, 708, 710, 712, 714, 716, 718, 720, 722, 724, 726, 728, 730, 732, 734, 736, 738, 740, 742, 744, 746, 748, 750, 752, 754, 756, 758, 760, 762, 764, 766, 768, 770, 772, 774, 776, 778, 780, 782, 784, 786, 788, 790, 792, 794, 796, 798, 800, 802, 804, 806, 808, 810, 812, 814, 816, 818, 820, 822, 824, 826, 828, 830, 832, 834, 836, 838, 840, 842, 844, 846, 848, 850, 852, 854, 856, 858, 860, 862, 864, 866, 868, 870, 872, 874, 876, 878, 880, 882, 884, 886, 888, 890, 892, 894, 896, 898, 900, 902, 904, 906, 908, 910, 912, 914, 916, 918, 920, 922, 924, 926, 928, 930, 932, 934, 936, 938, 940, 942, 944, 946, 948, 950, 952, 954, 956, 958, 960, 962, 964, 966, 968, 970, 972, 974, 976, 978, 980, 982, 984, 986, 988, 990, 992, 994, 996, 998, 1000
    
    # Array size
    spectrum_size: .long 150
    
    # Graph adjacency list structure
    graph_nodes: .long 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    graph_edges: .long 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    graph_count: .long 0
    
    # Buffer for output formatting
    buffer: .space 256

.section .text
.globl _start

_start:
    # Initialize registers
    movl spectrum_size, %ecx        # Load array size into ecx
    movl $0, %esi                   # Initialize index counter
    movl $0, %edi                   # Initialize edge counter
    
    # Outer loop: for each mass in spectrum
outer_loop:
    cmpl $0, %ecx                   # Check if we've processed all elements
    jz end_program                  # If yes, exit
    
    # Load current mass
    movl spectrum(,%esi,4), %eax    # Load current spectrum value into eax
    
    # Inner loop: compare with all other masses
    movl $0, %edx                   # Reset inner counter
    
inner_loop:
    cmpl %ecx, %edx                 # Check if we've compared with all elements
    jge next_outer                  # If yes, move to next outer element
    
    # Skip comparing with itself
    cmpl %esi, %edx
    je skip_inner
    
    # Load second mass
    movl spectrum(,%edx,4), %ebx    # Load second spectrum value into ebx
    
    # Calculate difference (mass1 - mass2)
    subl %ebx, %eax                 # Calculate difference
    
    # Store the edge in graph structure
    pushl %eax                      # Push difference (weight)
    pushl %esi                      # Push source index
    pushl %edx                      # Push destination index
    
    # Store in adjacency list format
    movl %eax, graph_edges(,%edi,4) # Store weight
    addl $1, %edi                   # Increment edge counter
    
skip_inner:
    incl %edx                       # Increment inner counter
    jmp inner_loop                  # Continue inner loop

next_outer:
    incl %esi                       # Increment outer counter
    decl %ecx                       # Decrement array size counter
    jmp outer_loop                  # Continue outer loop

end_program:
    # Exit program
    movl $1, %eax                   # sys_exit
    movl $0, %ebx                   # exit status
    int $0x80                       # system call
```

## Explanation of Key Components

### Data Section
- **spectrum**: Array containing mass values from the spectrum
- **spectrum_size**: Number of elements in the spectrum array
- **graph_nodes/edges**: Structures to store graph information
- **buffer**: Space for output formatting

### Algorithm Steps
1. **Outer Loop**: Iterate through each mass in the spectrum
2. **Inner Loop**: For each mass, compare with all other masses
3. **Difference Calculation**: Compute mass differences (source - destination)
4. **Graph Construction**: Store edges in adjacency list format

### Key Features
- Uses nested loops to generate all possible transitions
- Stores edge weights (differences) in the graph structure
- Handles all pairs of masses exactly once
- Produces a directed graph representation suitable for further processing

This assembly implementation constructs the graph where each directed edge represents a mass transition with weight equal to the difference between the two masses.