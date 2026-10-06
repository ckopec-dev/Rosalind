# Rosalind Problem: Implement DistanceBetweenPatternAndStrings

## Problem Understanding

The task is to find the distance between a pattern (k-mer) and a collection of DNA strings. The distance is defined as the sum of the minimum Hamming distances from the pattern to each string in the collection.

## Assembly Implementation

```assembly
.data
    # Input pattern and DNA strings
    pattern: .ascii "ACGT"      # Example pattern (4-mer)
    pattern_len: .long 4
    
    dna_strings: .ascii "AAATTTCCGGGTTT"
    dna_strings_len: .long 18
    
    # Arrays for storing string lengths
    string_lengths: .long 6, 6, 6   # Lengths of three strings
    num_strings: .long 3            # Number of DNA strings
    
    # Result storage
    min_distance: .long 0
    total_distance: .long 0
    
    # Temporary variables
    temp_dist: .long 0
    min_hamming: .long 1000000
    
.text
.globl _start

_start:
    # Initialize registers
    movl pattern_len, %ecx          # Pattern length
    movl num_strings, %edx          # Number of strings
    
    # Main loop through each DNA string
    xorl %esi, %esi                 # String index counter (i = 0)
    
string_loop:
    # Check if we've processed all strings
    cmpl %edx, %esi
    jge end_program
    
    # Get current string length
    movl string_lengths(,%esi,4), %eax
    
    # Find minimum Hamming distance between pattern and current string
    call find_min_hamming_distance
    
    # Add to total distance
    addl temp_dist, total_distance
    
    # Increment string counter
    incl %esi
    jmp string_loop

# Function: find_min_hamming_distance
# Input: pattern in memory, current_string in memory, pattern_len in %ecx
# Output: minimum Hamming distance in temp_dist
find_min_hamming_distance:
    pushl %ebp
    movl %esp, %ebp
    
    # Initialize variables
    xorl %edi, %edi                 # Position in string (j = 0)
    movl $1000000, min_hamming      # Initialize to large number
    
string_position_loop:
    # Check if we've gone beyond possible alignments
    cmpl %eax, %edi
    jge calculate_min_distance
    
    # Calculate Hamming distance for current alignment
    call calculate_hamming_distance
    
    # Update minimum if current is smaller
    cmpl temp_dist, min_hamming
    jge skip_update
    movl temp_dist, min_hamming
    
skip_update:
    incl %edi
    jmp string_position_loop

calculate_min_distance:
    movl min_hamming, temp_dist
    popl %ebp
    ret

# Function: calculate_hamming_distance
# Input: pattern in memory, current_string in memory, pattern_len in %ecx, position in %edi
# Output: Hamming distance in temp_dist
calculate_hamming_distance:
    pushl %ebp
    movl %esp, %ebp
    
    xorl %esi, %esi                 # Character index (k = 0)
    xorl %eax, %eax                 # Distance counter
    
hamming_char_loop:
    cmpl %ecx, %esi                 # Check if we've compared all characters
    jge hamming_done
    
    # Get pattern character at position k
    movb pattern(,%esi,1), %bl
    
    # Get string character at current alignment + k
    movb current_string(%edi,%esi,1), %cl
    
    # Compare characters
    cmpb %cl, %bl
    je skip_increment               # If equal, don't increment distance
    
    incl %eax                       # Increment Hamming distance
    
skip_increment:
    incl %esi
    jmp hamming_char_loop

hamming_done:
    movl %eax, temp_dist
    popl %ebp
    ret

# Function: get_string_length
# Input: pointer to string in %eax
# Output: length in %ebx
get_string_length:
    pushl %ebp
    movl %esp, %ebp
    
    xorl %ebx, %ebx                 # Initialize counter
    
length_loop:
    movb (%eax,%ebx,1), %cl
    cmpb $0, %cl                    # Check for null terminator
    je length_done
    
    incl %ebx
    jmp length_loop
    
length_done:
    movl %ebx, %eax
    popl %ebp
    ret

end_program:
    # Store final result
    movl total_distance, min_distance
    
    # Exit program
    movl $1, %eax                   # sys_exit
    movl $0, %ebx                   # exit status
    int $0x80
```

## Alternative Implementation (More Efficient)

```assembly
.data
    # Pattern and DNA strings
    pattern: .ascii "ACGT"
    pattern_len: .long 4
    
    # Three DNA strings of equal length for this example
    string1: .ascii "AAATTT"
    string2: .ascii "CCGGTT" 
    string3: .ascii "AATCGG"
    
    # String lengths (assuming all same length)
    str_len: .long 6
    
    # Results
    total_dist: .long 0

.text
.globl _start

_start:
    movl $pattern, %esi             # Load pattern address
    movl pattern_len, %ecx          # Pattern length
    
    # Process first string
    movl $string1, %edi
    call hamming_distance
    addl %eax, total_dist
    
    # Process second string  
    movl $string2, %edi
    call hamming_distance
    addl %eax, total_dist
    
    # Process third string
    movl $string3, %edi
    call hamming_distance
    addl %eax, total_dist
    
    # Exit
    movl $1, %eax
    movl $0, %ebx
    int $0x80

# Function: hamming_distance
# Input: pattern in %esi, string in %edi, length in %ecx
# Output: Hamming distance in %eax
hamming_distance:
    xorl %eax, %eax                 # Distance counter = 0
    xorl %edx, %edx                 # Character index
    
distance_loop:
    cmpl %ecx, %edx                 # Check if done with pattern
    jge distance_done
    
    # Compare characters at current position
    movb (%esi,%edx,1), %bl         # Load pattern char
    movb (%edi,%edx,1), %cl         # Load string char
    
    cmpb %cl, %bl                   # Compare characters
    jne increment_distance          # If different, increment distance
    
    jmp next_char                   # Continue to next character
    
increment_distance:
    incl %eax                       # Increment Hamming distance
    
next_char:
    incl %edx                       # Move to next character
    jmp distance_loop
    
distance_done:
    ret
```

## Key Implementation Details

1. **Pattern-String Comparison**: The algorithm compares the pattern with all possible substrings of each DNA string
2. **Hamming Distance Calculation**: For each alignment, count mismatches between characters
3. **Minimum Finding**: Find the minimum Hamming distance for each string
4. **Summation**: Add up minimum distances from all strings

## Time Complexity
- O(n × m × k) where n = number of strings, m = length of each string, k = pattern length

## Space Complexity  
- O(1) additional space (excluding input storage)

This assembly implementation handles the core algorithmic requirements of finding the distance between a pattern and multiple DNA strings using Hamming distance calculations.