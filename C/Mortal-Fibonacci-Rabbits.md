# Rosalind Problem: Mortal Fibonacci Rabbits

## Problem Description
In nature, rabbits live for a fixed number of months before dying. We are given the number of months `n` and the lifespan of rabbits `m`. We start with one pair of newborn rabbits, and each pair produces a new pair every month starting from the second month. Each pair of rabbits lives for exactly `m` months (including the month they were born). How many pairs of rabbits will there be after `n` months?

## Solution Approach
This is an extension of the classic Fibonacci problem where we need to account for rabbit mortality. We'll use dynamic programming with a sliding window approach, keeping track of rabbits in different age groups.

## C Implementation

```c
#include <stdio.h>
#include <stdlib.h>

long long mortal_fibonacci_rabbits(int n, int m) {
    // Edge cases
    if (n == 0) return 0;
    if (n == 1 || n == 2) return 1;
    
    // Array to store the number of rabbits in each age group
    // Index 0 represents rabbits that just reproduced (age 1)
    // Index m-1 represents rabbits that are about to die (age m)
    long long *rabbits = (long long*)calloc(m, sizeof(long long));
    
    // Initially we have one pair of newborn rabbits
    rabbits[0] = 1;
    
    // Simulate each month
    for (int i = 2; i <= n; i++) {
        long long new_pairs = 0;
        
        // All rabbits except the oldest ones reproduce
        for (int j = 1; j < m; j++) {
            new_pairs += rabbits[j];
        }
        
        // Shift all age groups forward by one month
        for (int j = m - 1; j > 0; j--) {
            rabbits[j] = rabbits[j - 1];
        }
        
        // New pairs are born and placed in the youngest age group
        rabbits[0] = new_pairs;
    }
    
    // Sum all rabbits from all age groups
    long long total = 0;
    for (int i = 0; i < m; i++) {
        total += rabbits[i];
    }
    
    free(rabbits);
    return total;
}

int main() {
    int n, m;
    
    // Read input
    printf("Enter number of months (n): ");
    scanf("%d", &n);
    printf("Enter lifespan of rabbits (m): ");
    scanf("%d", &m);
    
    // Solve and print result
    long long result = mortal_fibonacci_rabbits(n, m);
    printf("Number of rabbit pairs after %d months: %lld\n", n, result);
    
    return 0;
}
```

## How It Works

1. **Initialization**: Create an array where each index represents rabbits in a specific age group (0 to m-1)
2. **Base Case**: Start with one pair of newborn rabbits (age 0)
3. **Monthly Simulation**: 
   - Calculate new pairs born from all reproductive rabbits (ages 1 to m-2)
   - Shift all rabbits forward by one month (aging)
   - Place new pairs in the youngest age group
4. **Result**: Sum all rabbits across all age groups

## Example Walkthrough
For n=6, m=3:
- Month 1: 1 pair (age 0)
- Month 2: 1 pair (age 0) 
- Month 3: 1 pair (age 0) + 1 pair (from previous generation) = 2 pairs
- Month 4: New pairs from age 1 + age 2 rabbits = 1 pair + 1 pair = 2 pairs

## Time and Space Complexity
- **Time Complexity**: O(n × m)
- **Space Complexity**: O(m)

## Sample Input/Output
```
Input:
n = 6, m = 3

Output:
Number of rabbit pairs after 6 months: 4
```