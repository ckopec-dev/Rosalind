# Rosalind Problem: Rabbits and Recurrence Relations (C Solution)

## Problem Understanding

This is the classic Fibonacci rabbit population problem where:
- Each pair of rabbits produces a new pair every month
- Rabbits mature after 1 month and start reproducing
- We need to calculate the number of rabbit pairs after n months with k offspring per pair

## Solution

```c
#include <stdio.h>

long long fibonacci_rabbits(int n, int k) {
    // Base cases
    if (n == 1 || n == 2) {
        return 1;
    }
    
    // For n > 2, F(n) = F(n-1) + k*F(n-2)
    long long prev2 = 1;  // F(1)
    long long prev1 = 1;  // F(2)
    long long current;
    
    for (int i = 3; i <= n; i++) {
        current = prev1 + k * prev2;
        prev2 = prev1;
        prev1 = current;
    }
    
    return current;
}

int main() {
    int n, k;
    
    // Read input values
    printf("Enter number of months (n): ");
    scanf("%d", &n);
    printf("Enter number of offspring per pair (k): ");
    scanf("%d", &k);
    
    // Calculate and print result
    long long result = fibonacci_rabbits(n, k);
    printf("Number of rabbit pairs after %d months: %lld\n", n, result);
    
    return 0;
}
```

## Alternative Implementation with Array Storage

```c
#include <stdio.h>

long long fibonacci_rabbits_array(int n, int k) {
    if (n <= 0) return 0;
    if (n == 1 || n == 2) return 1;
    
    long long fib[n + 1];
    fib[1] = 1;
    fib[2] = 1;
    
    for (int i = 3; i <= n; i++) {
        fib[i] = fib[i-1] + k * fib[i-2];
    }
    
    return fib[n];
}

int main() {
    int n, k;
    
    printf("Enter number of months (n): ");
    scanf("%d", &n);
    printf("Enter number of offspring per pair (k): ");
    scanf("%d", &k);
    
    long long result = fibonacci_rabbits_array(n, k);
    printf("Number of rabbit pairs after %d months: %lld\n", n, result);
    
    return 0;
}
```

## Sample Input/Output

**Input:**
```
Enter number of months (n): 5
Enter number of offspring per pair (k): 3
```

**Output:**
```
Number of rabbit pairs after 5 months: 19
```

## Explanation

The recurrence relation is:
- F(1) = 1, F(2) = 1
- F(n) = F(n-1) + k × F(n-2) for n > 2

Where:
- n = number of months
- k = number of offspring produced by each pair per month
- F(n) = number of rabbit pairs in month n

For the example (n=5, k=3):
- Month 1: 1 pair
- Month 2: 1 pair  
- Month 3: 1 + 3×1 = 4 pairs
- Month 4: 4 + 3×1 = 7 pairs
- Month 5: 7 + 3×4 = 19 pairs

## Time and Space Complexity

- **Time Complexity**: O(n)
- **Space Complexity**: O(1) for iterative approach, O(n) for array approach

The iterative solution is more memory-efficient and recommended for large inputs.