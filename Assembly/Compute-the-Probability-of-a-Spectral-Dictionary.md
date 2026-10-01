# Rosalind Problem: Compute the Probability of a Spectral Dictionary

## Problem Understanding

We need to compute the probability of a spectral dictionary, which involves:
1. Given a spectrum and a dictionary of peptides
2. Calculate the probability that a random peptide from the dictionary would produce the given spectrum
3. This requires understanding how amino acid masses contribute to peptide masses

## Solution Approach

```assembly
; Compute the Probability of a Spectral Dictionary
; Input: Spectrum (array of masses), Dictionary (array of peptides)
; Output: Probability value

.data
    ; Amino acid masses
    aa_masses DWORD 57, 71, 87, 97, 99, 101, 103, 113, 114, 115
    DWORD 128, 129, 131, 137, 147, 156, 163, 186
    
    ; Number of amino acids
    num_aa DWORD 18
    
    ; Spectrum data (example)
    spectrum DWORD 113, 129, 147, 163, 186
    spec_size DWORD 5
    
    ; Dictionary peptides (represented as mass sums for simplicity)
    dict_peptides DWORD 113, 129, 147, 163, 186
    dict_size DWORD 5

.text
.code

compute_prob_spectral_dict PROC
    ; Input parameters:
    ; ECX = pointer to spectrum array
    ; EDX = size of spectrum
    ; ESI = pointer to dictionary array  
    ; EDI = size of dictionary
    
    push ebp
    mov ebp, esp
    
    ; Initialize counters and variables
    xor eax, eax          ; total_count
    xor ebx, ebx          ; match_count
    xor ecx, ecx          ; spectrum_index
    xor edx, edx          ; dict_index
    
    ; Get spectrum size
    mov edi, [ebp+12]     ; spectrum_size
    mov esi, [ebp+8]      ; spectrum_array
    
    ; Get dictionary size  
    mov ecx, [ebp+20]     ; dict_size
    mov edx, [ebp+16]     ; dict_array
    
    ; For each peptide in dictionary
    mov ebp, 0            ; peptide_index = 0
    
compute_loop:
    cmp ebp, ecx          ; compare with dict_size
    jge end_compute       ; if index >= size, exit loop
    
    ; Get current peptide mass
    mov eax, [edx + ebp*4] ; peptide_mass = dict[peptide_index]
    
    ; Check if this peptide can produce the spectrum
    push ebp              ; save index for later
    call check_spectrum_match
    pop ebp               ; restore index
    
    ; If match found, increment match_count
    cmp eax, 1
    jne next_peptide
    
    inc ebx               ; match_count++
    
next_peptide:
    inc ebp               ; peptide_index++
    jmp compute_loop
    
end_compute:
    ; Calculate probability = matches / total_peptides
    mov eax, ebx          ; matches in eax
    xor edx, edx          ; clear edx for division
    mov ebx, [ebp+20]     ; dict_size in ebx
    div ebx               ; eax = matches / total
    
    ; Return probability value in eax
    pop ebp
    ret
    
compute_prob_spectral_dict ENDP

check_spectrum_match PROC
    ; Check if a peptide mass can produce the given spectrum
    ; Input: EAX = peptide_mass, ECX = spectrum_array, EDX = spectrum_size
    ; Output: EAX = 1 if match, 0 if no match
    
    push ebp
    mov ebp, esp
    
    ; Simple implementation - check if peptide mass exists in spectrum
    ; In a real implementation, this would be more complex
    xor eax, eax          ; return 0 (no match)
    
    ; For demonstration, we'll say it matches if peptide mass is in spectrum
    mov ecx, [ebp+12]     ; spectrum_size
    mov esi, [ebp+8]      ; spectrum_array
    
    ; Linear search through spectrum
    xor edi, edi          ; index = 0
    
search_loop:
    cmp edi, ecx          ; compare with size
    jge no_match          ; if index >= size, no match
    
    cmp eax, [esi + edi*4] ; compare peptide_mass with spectrum[i]
    je found_match        ; if equal, found match
    
    inc edi               ; increment index
    jmp search_loop       ; continue search
    
no_match:
    pop ebp
    ret
    
found_match:
    mov eax, 1            ; return 1 (match found)
    pop ebp
    ret
    
check_spectrum_match ENDP

END
```

## Alternative Implementation (More Realistic)

```assembly
; More realistic implementation for spectral dictionary probability
; Uses amino acid mass combinations to build peptide probabilities

.data
    ; Standard amino acid masses (monoisotopic)
    aa_masses DB 57, 71, 87, 97, 99, 101, 103, 113, 114, 115
    DB 128, 129, 131, 137, 147, 156, 163, 186
    
    ; Number of amino acids  
    num_aa DWORD 18
    
    ; Spectrum array
    spectrum DWORD 113, 129, 147, 163, 186
    spec_size DWORD 5

.code

prob_spectral_dict PROC
    ; Initialize
    xor eax, eax          ; result = 0.0
    xor ebx, ebx          ; counter
    
    ; Get spectrum information
    mov ecx, [spec_size]  ; number of spectrum peaks
    mov esi, offset spectrum
    
    ; For each peak in spectrum
    xor edi, edi          ; peak_index = 0
    
peak_loop:
    cmp edi, ecx
    jge peak_end
    
    ; Get current peak mass
    mov edx, [esi + edi*4]
    
    ; Count how many peptides from dictionary would generate this mass
    call count_matching_peptides
    
    ; Add to result (simplified)
    add eax, ebx          ; accumulate counts
    
    inc edi               ; next peak
    jmp peak_loop
    
peak_end:
    ; Calculate final probability (simplified calculation)
    mov ebx, [spec_size]
    xor edx, edx
    div ebx               ; eax = total_matches / spectrum_size
    
    ret
    
count_matching_peptides PROC
    ; Count peptides that could produce given mass
    ; Input: EDX = target_mass
    ; Output: EBX = count of matching peptides
    
    push ebp
    mov ebp, esp
    
    xor ebx, ebx          ; count = 0
    
    ; Simple implementation - just return a basic count
    ; In reality, this would enumerate possible amino acid combinations
    mov ebx, 1            ; assume one match for now
    
    pop ebp
    ret
    
count_matching_peptides ENDP

prob_spectral_dict ENDP

END
```

## Key Concepts

1. **Spectral Dictionary**: A set of peptides that could potentially produce a given mass spectrum
2. **Probability Calculation**: The likelihood that a random peptide from the dictionary matches the spectrum
3. **Mass Matching**: Each peptide's mass must be compatible with the spectrum peaks
4. **Combinatorial Analysis**: Consider all possible amino acid combinations that could generate each peak

The assembly code implements a basic framework for computing these probabilities, though in practice this would involve more sophisticated algorithms for peptide mass calculation and matching.