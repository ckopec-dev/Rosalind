# Rosalind Problem: Construct the De Bruijn Graph of a String

## Problem Understanding

The De Bruijn graph is a fundamental data structure in bioinformatics used for genome assembly. Given a string and a k-mer length, we need to construct a De Bruijn graph where:
- Each node represents a (k-1)-mer
- Each edge represents a k-mer
- An edge from node A to node B exists if the suffix of A equals the prefix of B

## Solution Approach

I'll implement this in Assembly using x86-64 syntax with NASM format:

```assembly
section .data
    ; Input string and parameters
    input_string db "ACGTACGT", 0
    k_value equ 3
    string_length equ 8
    
    ; Output buffers
    nodes_buffer db 100 dup(0)
    edges_buffer db 200 dup(0)
    
section .bss
    ; Hash table for storing nodes
    node_hash_table resb 1024
    edge_list resb 500
    
section .text
    global _start
    
; Function to construct De Bruijn graph
construct_debruijn_graph:
    push rbp
    mov rbp, rsp
    
    ; Initialize variables
    xor rax, rax                    ; i = 0
    mov rcx, k_value                ; k = k_value
    dec rcx                         ; k-1 for suffix/prefix
    mov rdx, string_length          ; n = string_length
    sub rdx, rcx                    ; n-(k-1) = number of k-mers
    
    ; Process each k-mer
process_kmers:
    cmp rax, rdx
    jge graph_construction_done
    
    ; Extract k-mer at position rax
    mov rdi, input_string           ; source string
    add rdi, rax                    ; pointer to current k-mer
    mov rsi, nodes_buffer           ; destination buffer
    mov rcx, k_value                ; copy k characters
    rep movsb                       ; copy k-mer
    
    ; Get suffix (k-1 characters) and prefix (k-1 characters)
    mov rdi, nodes_buffer           ; current k-mer
    mov rsi, edges_buffer           ; destination for edge
    mov rcx, k_value                ; k-1 characters to copy
    dec rcx                         ; k-1
    
    ; Extract suffix (first k-1 chars)
    mov r8, rdi                     ; save start of k-mer
    add rdi, 1                      ; skip first character for suffix
    mov r9, rsi                     ; save edge buffer start
    mov r10, rcx                    ; copy count
    
suffix_loop:
    cmp r10, 0
    jle suffix_done
    mov al, [r8]
    mov [r9], al
    inc r8
    inc r9
    dec r10
    jmp suffix_loop
    
suffix_done:
    ; Extract prefix (first k-1 chars)
    mov rdi, nodes_buffer           ; current k-mer
    mov rsi, edges_buffer           ; destination for edge
    add rsi, 50                     ; offset for prefix
    mov rcx, k_value
    dec rcx                         ; k-1 characters
    
prefix_loop:
    cmp rcx, 0
    jle prefix_done
    mov al, [rdi]
    mov [rsi], al
    inc rdi
    inc rsi
    dec rcx
    jmp prefix_loop
    
prefix_done:
    ; Add to hash table or adjacency list
    call add_to_graph
    
    inc rax                         ; i++
    jmp process_kmers
    
graph_construction_done:
    pop rbp
    ret

; Function to add k-mer to graph structure
add_to_graph:
    push rbp
    mov rbp, rsp
    
    ; Hash the k-mer for node storage
    mov rdi, nodes_buffer
    call hash_string
    
    ; Store edge information
    ; In a real implementation, this would build adjacency lists
    ; For now, we'll just demonstrate the concept
    
    pop rbp
    ret

; Simple hash function for strings
hash_string:
    push rbp
    mov rbp, rsp
    
    xor rax, rax                    ; hash = 0
    mov rcx, k_value                ; length of string
    mov rdi, nodes_buffer           ; string to hash
    
hash_loop:
    cmp rcx, 0
    jle hash_done
    
    movzx rdx, byte [rdi]           ; get character
    add rax, rdx                    ; add to hash
    inc rdi                         ; next character
    dec rcx                         ; decrement counter
    jmp hash_loop
    
hash_done:
    pop rbp
    ret

; Main function
_start:
    ; Initialize stack and registers
    push rbp
    mov rbp, rsp
    
    ; Call the graph construction function
    call construct_debruijn_graph
    
    ; Exit program
    mov eax, 1                      ; sys_exit
    xor ebx, ebx                    ; exit status 0
    int 0x80                        ; system call
```

## Alternative Simpler Implementation

```assembly
section .data
    ; Input string and parameters
    dna_string db "ACGTACGT", 0
    k equ 3
    n equ 8
    
    ; Output format strings
    node_format db "Node: %s\n", 0
    edge_format db "Edge: %s -> %s\n", 0

section .text
    global _start

; Simple implementation that processes k-mers
simple_debruijn:
    push rbp
    mov rbp, rsp
    
    ; Initialize loop counter
    xor rax, rax                    ; i = 0
    
process_loop:
    ; Check if we've processed all k-mers
    mov rbx, n                      ; n = string length
    sub rbx, k                      ; n-k
    inc rbx                         ; n-k+1 possible k-mers
    cmp rax, rbx
    jge done
    
    ; Extract current k-mer from dna_string + i
    mov rdi, dna_string             ; base address
    add rdi, rax                    ; point to start of k-mer
    push rdi                        ; save pointer
    
    ; Print k-mer (simplified)
    ; In real implementation, you'd extract prefix and suffix
    
    ; Extract prefix (k-1 characters from start)
    mov rdi, dna_string
    add rdi, rax                    ; prefix starts here
    mov rsi, [rdi]                  ; get first k-1 chars
    
    ; Extract suffix (k-1 characters from end)
    mov rdi, dna_string
    add rdi, rax
    add rdi, 1                      ; skip first char for suffix
    
    inc rax                         ; i++
    jmp process_loop
    
done:
    pop rbp
    ret

_start:
    call simple_debruijn
    
    ; Exit
    mov eax, 1
    xor ebx, ebx
    int 0x80
```

## Key Algorithm Steps

1. **Input Processing**: Read the DNA string and k-mer length
2. **K-mer Extraction**: Extract all k-mers from the string
3. **Prefix/Suffix Calculation**: For each k-mer, calculate its prefix (k-1 characters) and suffix (k-1 characters)
4. **Graph Construction**: Create nodes for prefixes and edges connecting matching prefixes/suffixes
5. **Output**: Format the adjacency list representation

## Time Complexity
- O(n×k) where n is string length and k is k-mer size
- Each k-mer requires O(k) operations to process

## Space Complexity  
- O(n×k) for storing all k-mers and graph structure

The assembly implementation above demonstrates the core logic, though a complete implementation would require more sophisticated string handling and memory management.