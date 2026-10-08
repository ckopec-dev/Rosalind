# Rosalind Problem: Reconstruct a String from its Genome Path

## Problem Statement
Given a sequence of k-mers that form a genome path, reconstruct the original DNA string. The genome path consists of a series of overlapping k-mers where each k-mer overlaps with the next one by k-1 nucleotides.

## Solution in Assembly

```assembly
; Reconstruct a string from its genome path
; Input: List of k-mers (strings of length k)
; Output: Original DNA string reconstructed from the path

.data
    ; Input k-mers (example)
    kmers db "ACGT", "CGTA", "GTAC", "TACG", 0  ; Null-terminated list
    max_kmers equ 100
    k_value equ 4
    
    ; Output buffer
    output_string db 200 dup(?)
    string_length dd 0

.code
main proc
    ; Initialize registers
    mov esi, offset kmers    ; Pointer to first k-mer
    mov edi, offset output_string  ; Pointer to output buffer
    xor ecx, ecx             ; Counter for k-mers
    
    ; Get first k-mer and copy it to output
    mov eax, esi
    mov ebx, 0
    mov byte ptr [edi], byte ptr [eax]  ; Copy first character
    inc edi
    
    ; Process remaining k-mers
    process_loop:
        ; Check if we've processed all k-mers
        cmp ecx, max_kmers
        jge end_process
        
        ; Get current k-mer
        mov eax, esi
        add eax, k_value        ; Move to next k-mer
        mov esi, eax
        
        ; Copy k-1 characters from the beginning of next k-mer
        ; (the first character is already covered by overlap)
        mov eax, esi
        mov ebx, 1              ; Start from second character (index 1)
        mov edx, k_value
        
        copy_loop:
            cmp ebx, edx
            jge copy_done
            
            mov al, byte ptr [eax + ebx]  ; Get character at position ebx
            mov byte ptr [edi], al
            inc edi
            inc ebx
            jmp copy_loop
            
        copy_done:
            inc ecx
            jmp process_loop
    
    end_process:
        ; Null terminate the output string
        mov byte ptr [edi], 0
        
        ; Return result (output_string)
        ret
main endp

; Alternative implementation using a more efficient approach
reconstruct_genome proc
    ; Input: pointer to array of k-mers
    ; Output: reconstructed DNA string
    
    push ebp
    mov ebp, esp
    
    ; Get parameters
    mov eax, [ebp + 8]      ; Pointer to kmers array
    mov ebx, [ebp + 12]     ; Number of k-mers
    
    ; Initialize output buffer
    xor ecx, ecx            ; Counter
    mov edi, offset output_string
    
    ; Copy first k-mer completely
    mov esi, eax
    mov edx, k_value        ; Length of each k-mer
    copy_first:
        cmp ecx, edx
        jge first_done
        
        mov al, byte ptr [esi + ecx]
        mov byte ptr [edi + ecx], al
        inc ecx
        jmp copy_first
    
    first_done:
        ; Process remaining k-mers (only last k-1 characters)
        mov ecx, 1              ; Start from second k-mer
        mov esi, eax            ; Reset pointer
        
    process_remaining:
        cmp ecx, ebx
        jge done_processing
        
        ; Move to next k-mer in array
        add esi, k_value        ; Skip to next k-mer
        
        ; Copy only last (k-1) characters
        mov edx, 1              ; Start from index 1
        mov edi, offset output_string
        add edi, ecx            ; Position in output string
        add edi, k_value        ; Move to correct position
        
    copy_last_chars:
        cmp edx, k_value
        jge last_done
        
        mov al, byte ptr [esi + edx]
        mov byte ptr [edi], al
        inc edi
        inc edx
        jmp copy_last_chars
    
    last_done:
        inc ecx
        jmp process_remaining
    
    done_processing:
        ; Null terminate
        mov byte ptr [edi], 0
        mov eax, offset output_string
        
    pop ebp
    ret
reconstruct_genome endp

end
```

## Algorithm Explanation

1. **Input Processing**: Read the sequence of k-mers from input
2. **First K-mer**: Copy the entire first k-mer to the result string
3. **Subsequent K-mers**: For each remaining k-mer, copy only the last (k-1) characters 
4. **Overlap Handling**: Since consecutive k-mers overlap by (k-1) nucleotides, we only need to add the last character of each k-mer
5. **Output**: Constructed DNA string with proper overlaps

## Example

Input: `["ACGT", "CGTA", "GTAC", "TACG"]`
Output: `"ACGTACGT"`

The algorithm ensures that:
- First k-mer contributes all k characters
- Each subsequent k-mer contributes only 1 character (the last one) to maintain overlap
- Final string length = (number of k-mers × k) - (number of k-mers - 1)

## Time Complexity
O(n × k) where n is the number of k-mers and k is the k-mer length

## Space Complexity  
O(n × k) for storing the output string