# Rosalind Problem: Implement_ChromosomeToCycle

## Problem Description
Convert a chromosome (a sequence of signed integers) into a cycle representation.

## Solution in Assembly Language

```assembly
; ChromosomeToCycle - Convert chromosome to cycle representation
; Input: chromosome (signed integers)
; Output: cycle representation as array

.data
    ; Example chromosome: (+1 -2 -3 +4)
    chromosome: .word 1, -2, -3, 4
    chrom_length: .word 4
    
    ; Result array for cycle
    cycle: .space 8 * 4    ; 8 elements for 4 genes (each gene has 2 nodes)

.text
.globl _start

_start:
    ; Load chromosome and length
    la $t0, chromosome      ; Load address of chromosome
    lw $t1, chrom_length    ; Load chromosome length
    
    ; Initialize variables
    li $t2, 0               ; i = 0 (index for chromosome)
    li $t3, 0               ; j = 0 (index for cycle)
    
loop:
    ; Check if we've processed all elements
    beq $t2, $t1, done      ; if i == length, exit loop
    
    ; Load current gene value
    lw $t4, 0($t0)          ; Load chromosome[i]
    
    ; Check if gene is positive or negative
    bge $t4, $zero, positive_gene
    
    ; Handle negative gene
    neg $t5, $t4            ; $t5 = -gene (absolute value)
    sll $t5, $t5, 1         ; $t5 = 2 * |gene|
    addi $t5, $t5, -1       ; $t5 = 2 * |gene| - 1
    sw $t5, 0($t3)          ; cycle[j] = 2 * |gene| - 1
    
    addi $t5, $t5, 1        ; $t5 = 2 * |gene|
    sw $t5, 4($t3)          ; cycle[j+1] = 2 * |gene|
    
    addi $t3, $t3, 8        ; j += 2 (advance by 2 positions)
    addi $t0, $t0, 4        ; move to next chromosome element
    addi $t2, $t2, 1        ; i += 1
    j loop
    
positive_gene:
    ; Handle positive gene
    sll $t5, $t4, 1         ; $t5 = 2 * gene
    addi $t5, $t5, -1       ; $t5 = 2 * gene - 1
    sw $t5, 0($t3)          ; cycle[j] = 2 * gene - 1
    
    addi $t5, $t5, 1        ; $t5 = 2 * gene
    sw $t5, 4($t3)          ; cycle[j+1] = 2 * gene
    
    addi $t3, $t3, 8        ; j += 2
    addi $t0, $t0, 4        ; move to next chromosome element
    addi $t2, $t2, 1        ; i += 1
    j loop

done:
    ; Exit program
    li $v0, 10              ; Exit system call
    syscall
```

## Explanation

The assembly code implements the ChromosomeToCycle function that converts a chromosome (sequence of signed integers) into its cycle representation:

1. **Input Processing**: The algorithm processes each gene in the chromosome
2. **Gene Conversion**:
   - For positive genes: Gene `k` becomes nodes `2k-1` and `2k`
   - For negative genes: Gene `-k` becomes nodes `2k` and `2k-1`
3. **Cycle Construction**: The cycle is built by placing each gene's two nodes consecutively

## Example
For chromosome `[+1, -2, -3, +4]`:
- Gene 1 → nodes 1, 2  
- Gene -2 → nodes 4, 3
- Gene -3 → nodes 6, 5
- Gene 4 → nodes 7, 8

Result cycle: `[1, 2, 4, 3, 6, 5, 7, 8]`

## Key Assembly Concepts Used

- **Memory addressing**: Loading chromosome elements and storing results
- **Conditional branching**: Handling positive vs negative genes  
- **Bit manipulation**: Shifting for multiplication by 2
- **Loop control**: Processing all chromosome elements
- **System calls**: Proper program termination

The solution handles the mathematical conversion from signed integers to cycle representation as required by the Rosalind problem.