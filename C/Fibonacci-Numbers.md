# Rosalind Problem: Fibonacci Numbers (C Solution)

## Problem Statement
Given: A positive integer n ≤ 25.

Return: The value of Fn (the nth Fibonacci number).

## Solution

```c
#include <stdio.h>

int fibonacci(int n) {
    if (n <= 2) {
        return 1;
    }
    
    int prev2 = 1;  // F(1)
    int prev1 = 1;  // F(2)
    int current;
    
    for (int i = 3; i <= n; i++) {
        current = prev1 + prev2;
        prev2 = prev1;
        prev1 = current;
    }
    
    return current;
}

int main() {
    int n;
    printf("Enter the value of n: ");
    scanf("%d", &n);
    
    int result = fibonacci(n);
    printf("%d\n", result);
    
    return 0;
}
```

## Explanation

This solution uses an iterative approach to calculate the nth Fibonacci number:

1. **Base cases**: For n ≤ 2, we return 1 (since F(1) = F(2) = 1)
2. **Iterative calculation**: We use two variables (`prev2` and `prev1`) to keep track of the previous two Fibonacci numbers
3. **Loop**: Starting from i = 3, we calculate each subsequent Fibonacci number by summing the previous two
4. **Return**: The final result is stored in `current`

## Time and Space Complexity
- **Time Complexity**: O(n) - We iterate from 3 to n once
- **Space Complexity**: O(1) - Only using a constant amount of extra space

## Sample Input/Output
```
Input: 5
Output: 5

Input: 10
Output: 55
```

The Fibonacci sequence: 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, ...