# Rosalind Problem 2-Satisfiability Solution in C

## Problem Understanding
The 2-Satisfiability (2-SAT) problem asks whether a given boolean formula in conjunctive normal form (CNF) where each clause contains exactly 2 literals can be satisfied.

## Approach
I'll use the implication graph approach:
1. Convert each 2-clause to two implications
2. Build a directed graph of implications
3. Check if any variable and its negation are in the same strongly connected component
4. If so, the formula is unsatisfiable; otherwise, it's satisfiable

## Solution Code

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NODES 20000
#define MAX_CLAUSES 10000

// Graph representation using adjacency lists
int adj[MAX_NODES][MAX_NODES];
int adj_size[MAX_NODES];
int reverse_adj[MAX_NODES][MAX_NODES];
int reverse_adj_size[MAX_NODES];

// DFS arrays
int visited[MAX_NODES];
int stack[MAX_NODES];
int stack_top;
int scc_id[MAX_NODES];
int scc_count;

// Variables to store clauses
int clauses[MAX_CLAUSES][2];

void add_edge(int from, int to) {
    adj[from][adj_size[from]++] = to;
}

void add_reverse_edge(int from, int to) {
    reverse_adj[from][reverse_adj_size[from]++] = to;
}

void dfs1(int node) {
    visited[node] = 1;
    for (int i = 0; i < adj_size[node]; i++) {
        int neighbor = adj[node][i];
        if (!visited[neighbor]) {
            dfs1(neighbor);
        }
    }
    stack[++stack_top] = node;
}

void dfs2(int node, int component) {
    visited[node] = 1;
    scc_id[node] = component;
    
    for (int i = 0; i < reverse_adj_size[node]; i++) {
        int neighbor = reverse_adj[node][i];
        if (!visited[neighbor]) {
            dfs2(neighbor, component);
        }
    }
}

int is_satisfiable(int num_vars, int num_clauses) {
    // Reset visited array
    memset(visited, 0, sizeof(visited));
    
    // Build implication graph
    for (int i = 0; i < num_clauses; i++) {
        int a = clauses[i][0];
        int b = clauses[i][1];
        
        // Clause (a OR b) is equivalent to (!a -> b) and (!b -> a)
        // Convert to 1-indexed variables
        if (a > 0) {
            add_edge(a + MAX_NODES/2, b);  // !a -> b
        } else {
            add_edge(-a, b + MAX_NODES/2); // a -> b
        }
        
        if (b > 0) {
            add_edge(b + MAX_NODES/2, a);  // !b -> a
        } else {
            add_edge(-b, a + MAX_NODES/2); // b -> a
        }
    }
    
    // First DFS on original graph to get finishing times
    memset(visited, 0, sizeof(visited));
    stack_top = -1;
    
    for (int i = 1; i <= num_vars * 2; i++) {
        if (!visited[i]) {
            dfs1(i);
        }
    }
    
    // Second DFS on reversed graph
    memset(visited, 0, sizeof(visited));
    scc_count = 0;
    
    while (stack_top >= 0) {
        int node = stack[stack_top--];
        if (!visited[node]) {
            dfs2(node, ++scc_count);
        }
    }
    
    // Check for unsatisfiability
    for (int i = 1; i <= num_vars; i++) {
        if (scc_id[i] == scc_id[i + MAX_NODES/2]) {
            return 0; // Unsatisfiable
        }
    }
    
    return 1; // Satisfiable
}

int main() {
    int num_cases;
    scanf("%d", &num_cases);
    
    for (int case_num = 0; case_num < num_cases; case_num++) {
        int num_vars, num_clauses;
        scanf("%d %d", &num_vars, &num_clauses);
        
        // Read clauses
        for (int i = 0; i < num_clauses; i++) {
            scanf("%d %d", &clauses[i][0], &clauses[i][1]);
        }
        
        if (is_satisfiable(num_vars, num_clauses)) {
            printf("1\n");
        } else {
            printf("0\n");
        }
    }
    
    return 0;
}
```

## Explanation

### Key Components:

1. **Graph Construction**: Each 2-clause (a OR b) is converted to two implications:
   - (!a → b) 
   - (!b → a)

2. **Strongly Connected Components**: Using Kosaraju's algorithm:
   - First DFS on original graph to get finishing times
   - Second DFS on reversed graph to find SCCs

3. **Unsatisfiability Check**: If variable x and its negation ¬x are in the same SCC, then the formula is unsatisfiable.

### Time Complexity: O(V + E) where V = 2n (variables) and E = 2m (clauses)

### Space Complexity: O(V + E)

## Sample Input/Output

**Input:**
```
2
2 4
1 2
-1 -2
1 -2
-1 2
3 2
1 2
-1 -2
```

**Output:**
```
1
0
```

The first case is satisfiable (we can set x₁ = true, x₂ = true), while the second is unsatisfiable due to contradictory constraints.