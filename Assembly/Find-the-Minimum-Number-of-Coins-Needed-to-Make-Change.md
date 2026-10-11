# Find the Minimum Number of Coins Needed to Make Change - Assembly Solution

## Problem Understanding

We need to find the minimum number of coins required to make a given amount of money using dynamic programming approach.

## Assembly Implementation

```assembly
.data
    coins:      .word 1, 3, 5        # Available coin denominations
    num_coins:  .word 3              # Number of coin types
    amount:     .word 11             # Target amount to make change
    
.text
.globl _start

_start:
    # Load parameters
    la $a0, coins            # Load address of coins array
    lw $a1, num_coins        # Load number of coins
    lw $a2, amount           # Load target amount
    
    # Allocate memory for dp array (amount + 1 words)
    addi $t0, $a2, 1         # Calculate size: amount + 1
    sll $t0, $t0, 2          # Multiply by 4 (word size)
    jal malloc               # Allocate memory
    
    move $s0, $v0            # Store dp array address in s0
    
    # Initialize dp[0] = 0
    sw $zero, 0($s0)
    
    # Initialize all other values to MAX_INT (0x7FFFFFFF)
    li $t1, 1                # Counter i = 1
    li $t2, 0x7FFFFFFF       # MAX_INT value
    
initialize_loop:
    bge $t1, $a2, dp_init_done
    sll $t3, $t1, 2          # Calculate offset: i * 4
    add $t4, $s0, $t3        # Calculate address of dp[i]
    sw $t2, 0($t4)           # Set dp[i] = MAX_INT
    addi $t1, $t1, 1         # i++
    j initialize_loop
    
dp_init_done:
    # Main DP algorithm
    li $t1, 1                # i = 1 (start from 1)
outer_loop:
    bgt $t1, $a2, program_exit
    
    # For each coin denomination
    li $t2, 0                # j = 0 (coin index)
inner_loop:
    bge $t2, $a1, inner_loop_done
    
    # Get current coin value
    sll $t3, $t2, 2          # Calculate offset: j * 4
    add $t4, $a0, $t3        # Calculate address of coins[j]
    lw $t5, 0($t4)           # Load coin value
    
    # Check if coin value <= current amount i
    bgt $t5, $t1, continue_loop
    
    # Calculate dp[i - coin_value]
    sub $t6, $t1, $t5        # Calculate (i - coin_value)
    sll $t7, $t6, 2          # Calculate offset: (i - coin_value) * 4
    add $t8, $s0, $t7        # Calculate address of dp[i - coin_value]
    lw $t9, 0($t8)           # Load dp[i - coin_value]
    
    # Check if dp[i - coin_value] is not MAX_INT
    beq $t9, $zero, continue_loop
    
    # Add 1 to get new value
    addi $t9, $t9, 1         # dp[i - coin_value] + 1
    
    # Get current dp[i] value
    sll $t7, $t1, 2          # Calculate offset: i * 4
    add $t8, $s0, $t7        # Calculate address of dp[i]
    lw $t10, 0($t8)          # Load current dp[i]
    
    # Compare and update minimum
    ble $t9, $t10, continue_loop
    
    sw $t9, 0($t8)           # Update dp[i] = min(dp[i], dp[i - coin_value] + 1)
    
continue_loop:
    addi $t2, $t2, 1         # j++
    j inner_loop
    
inner_loop_done:
    addi $t1, $t1, 1         # i++
    j outer_loop

program_exit:
    # Return minimum coins needed (dp[amount])
    sll $t0, $a2, 2          # Calculate offset: amount * 4
    add $t1, $s0, $t0        # Calculate address of dp[amount]
    lw $v0, 0($t1)           # Load result
    
    # Exit program
    li $v0, 1                # syscall: exit
    syscall

# Simple malloc implementation for demonstration
malloc:
    # This is a simplified version - in practice would use system calls
    li $v0, 9                # syscall: sbrk
    move $a0, $zero          # Allocate size = 0 (for demo)
    syscall
    jr $ra
```

## Key Algorithm Steps

1. **Initialization**: Create a DP array where `dp[i]` represents minimum coins needed for amount `i`
2. **Base Case**: `dp[0] = 0` (no coins needed for amount 0)
3. **Fill DP Array**: For each amount from 1 to target:
   - Try each coin denomination
   - If coin value ≤ current amount, check if using this coin gives better solution
   - Update `dp[i]` with minimum value found

## Time and Space Complexity

- **Time Complexity**: O(amount × number_of_coins)
- **Space Complexity**: O(amount)

## Example Execution

For coins = [1, 3, 5], amount = 11:
- dp[0] = 0
- dp[1] = 1 (using coin 1)
- dp[2] = 2 (using two coins of 1)
- dp[3] = 1 (using coin 3)
- ...
- dp[11] = 3 (using coins: 5 + 5 + 1 or 3 + 3 + 3 + 1 + 1)

The minimum number of coins needed is **3**.