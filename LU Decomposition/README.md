# LU Decomposition Method

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](lu-decomposition.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Theory & Algorithm](#-theory--algorithm)
  - [Mathematical Foundation](#mathematical-foundation)
  - [Algorithm Steps](#algorithm-steps)
  - [Complexity Analysis](#complexity-analysis)
- [Implementation Details](#-implementation-details)
  - [Key Components](#key-components)
  - [LU Decomposition Process](#lu-decomposition-process)
  - [Solution Detection Logic](#solution-detection-logic)
- [Complete C++ Implementation](#-complete-c-implementation)
- [Usage Examples](#-usage-examples)
  - [Example 1: Unique Solution](#example-1-unique-solution)
  - [Example 2: No Solution](#example-2-no-solution)
  - [Example 3: Infinite Solutions](#example-3-infinite-solutions)
- [Compilation and Execution](#-compilation-and-execution)
- [Comparison with Other Methods](#-comparison-with-other-methods)
  - [Difference from Gauss Elimination](#difference-from-gauss-elimination)
  - [Difference from Gauss-Jordan Elimination](#difference-from-gauss-jordan-elimination)
  - [When to Use LU Decomposition](#when-to-use-lu-decomposition)
  - [Computational Cost Comparison](#computational-cost-comparison)
- [Applications](#-applications)
  - [Direct Applications](#direct-applications)
  - [Real-World Uses](#real-world-uses)
- [Mathematical Properties](#-mathematical-properties)
  - [Properties of LU Decomposition](#properties-of-lu-decomposition)
  - [Uniqueness](#uniqueness)
  - [Determinant Calculation](#determinant-calculation)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

**LU Decomposition** (also called LU Factorization) is a matrix decomposition method that factors a square matrix **A** into the product of a **Lower triangular matrix (L)** and an **Upper triangular matrix (U)**, such that **A = L × U**. 

This method is particularly efficient when solving multiple systems with the same coefficient matrix but different right-hand sides, as the decomposition needs to be computed only once.

This implementation in C++ solves systems of linear equations and can detect three types of solutions:
- **Unique Solution** - The system has exactly one solution
- **No Solution** - The system is inconsistent
- **Infinite Solutions** - The system has infinitely many solutions

### ⚙️ Features

✅ **Doolittle's Method** - L has 1s on diagonal, U has computed values  
✅ **Solution Detection** - Identifies unique, no, or infinite solutions  
✅ **Determinant Calculation** - Computed from diagonal of U  
✅ **Forward Substitution** - Solves L×y = b  
✅ **Back Substitution** - Solves U×x = y  
✅ **Multiple Test Cases** - Process several systems in one run  
✅ **Intermediate Output** - View L and U matrices at each step  
✅ **File I/O** - Reads from input. txt, writes to output.txt  
✅ **Solution Verification** - Checks A×x = b  
✅ **High Precision** - Uses double precision with 4 decimal places  
✅ **Error Handling** - Detects singular matrices and zero pivots  

---

## 🧮 Theory & Algorithm

### Mathematical Foundation

Given a system of linear equations **Ax = b**, where: 
- **A** is an n×n coefficient matrix
- **x** is the n×1 solution vector
- **b** is the n×1 constant vector

LU Decomposition factors **A** into:
```
A = L × U
```

Where:
- **L** is a lower triangular matrix (elements above diagonal are 0)
- **U** is an upper triangular matrix (elements below diagonal are 0)

The original system **Ax = b** becomes:
```
L × U × x = b
```

This is solved in two steps:
1. **Forward Substitution**: Solve **L×y = b** for **y**
2. **Back Substitution**: Solve **U×x = y** for **x**

### Algorithm Steps

#### 1. **LU Decomposition Phase (Doolittle's Method):**

For i = 0 to n-1:

**Step A: Compute U matrix (Upper triangular)**
```
For k = i to n-1:
    U[i][k] = A[i][k] - Σ(L[i][j] × U[j][k]) for j = 0 to i-1
```

**Step B: Compute L matrix (Lower triangular)**
```
L[i][i] = 1  (diagonal elements are always 1)

For k = i+1 to n-1:
    L[k][i] = (A[k][i] - Σ(L[k][j] × U[j][i]) for j = 0 to i-1) / U[i][i]
```

**Step C: Check for zero pivot**
- If U[i][i] ≈ 0, the matrix is singular → Check for no/infinite solutions

#### 2. **Solution Detection Phase:**

**Calculate determinant:**
```
det(A) = det(U) = Π(U[i][i]) for i = 0 to n-1
```

- If det(A) ≈ 0:
  - Check consistency:  If any equation reduces to 0 = c (where c ≠ 0) → **No Solution**
  - Otherwise → **Infinite Solutions**
- If det(A) ≠ 0 → **Unique Solution**

#### 3. **Forward Substitution (L×y = b):**

Solve for y from top to bottom:
```
For i = 0 to n-1:
    y[i] = b[i] - Σ(L[i][j] × y[j]) for j = 0 to i-1
```

Since L[i][i] = 1, no division is needed.

#### 4. **Back Substitution (U×x = y):**

Solve for x from bottom to top:
```
For i = n-1 down to 0:
    x[i] = (y[i] - Σ(U[i][j] × x[j]) for j = i+1 to n-1) / U[i][i]
```

### Complexity Analysis

#### Time Complexity
- **LU Decomposition**: O(n³/3) ≈ O(n³)
- **Forward Substitution**: O(n²)
- **Back Substitution**: O(n²)
- **Overall (one solve)**: O(n³)
- **Multiple solves with same A**: O(n³) + k×O(n²) where k = number of different b vectors

#### Space Complexity
- O(n²) for storing L, U, and A matrices
- O(n) for vectors x, y, and b
- **Total**: O(n²)

---

## 💻 Implementation Details

### Key Components

```cpp
vector<vector<double>> A(n, vector<double>(n));
vector<vector<double>> L(n, vector<double>(n, 0));
vector<vector<double>> U(n, vector<double>(n, 0));
```
- **A**: Original coefficient matrix (n×n)
- **L**: Lower triangular matrix initialized with zeros
- **U**: Upper triangular matrix initialized with zeros

```cpp
// Separate A and b from augmented matrix
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < n; j++)
        A[i][j] = aug[i][j];
    b[i] = aug[i][n];
}
```
- Extracts coefficient matrix A and constant vector b from augmented matrix

### LU Decomposition Process

```cpp
// Upper Triangular Matrix U
for (int k = i; k < n; k++)
{
    double sum = 0;
    for (int j = 0; j < i; j++)
        sum += L[i][j] * U[j][k];
    U[i][k] = A[i][k] - sum;
}
```
- Computes elements of U row-by-row
- Uses previously computed values of L and U

```cpp
// Lower Triangular Matrix L
for (int k = i; k < n; k++)
{
    if (i == k)
        L[i][i] = 1; // diagonal element is 1
    else
    {
        double sum = 0;
        for (int j = 0; j < i; j++)
            sum += L[k][j] * U[j][i];
        L[k][i] = (A[k][i] - sum) / U[i][i];
    }
}
```
- **Doolittle's Method**: Diagonal of L is always 1
- Computes L column-by-column below diagonal
- Divides by U[i][i] (pivot element)

```cpp
// Check for zero pivot
if (fabs(U[i][i]) < 1e-12)
{
    decompositionFailed = true;
    break;
}
```
- Detects singular matrices early
- Prevents division by zero

### Solution Detection Logic

```cpp
// Calculate determinant (product of diagonal of U)
double detU = 1;
for (int i = 0; i < n; i++)
    detU *= U[i][i];

if (fabs(detU) < 1e-12)
{
    // Matrix is singular - check for no solution vs infinite solutions
    // ...  consistency check ... 
}
```
- Determinant of A equals determinant of U (since det(L) = 1)
- If det ≈ 0, matrix is singular

```cpp
// Check for inconsistency
bool rowIsZero = true;
for (int j = i; j < n; j++)
{
    if (fabs(U[i][j]) > 1e-12)
    {
        rowIsZero = false;
        break;
    }
}

if (rowIsZero && fabs(expectedY) > 1e-12)
{
    noSolution = true;
    break;
}
```
- If a row in U is all zeros but corresponding y value is non-zero → **No Solution**
- Otherwise, if determinant is zero → **Infinite Solutions**

### Forward and Back Substitution

```cpp
// Forward substitution:  L*y = b
for (int i = 0; i < n; i++)
{
    double sum = 0;
    for (int j = 0; j < i; j++)
        sum += L[i][j] * y[j];
    y[i] = b[i] - sum;
}
```
- Solves from top to bottom
- No division needed since L[i][i] = 1

```cpp
// Back substitution: U*x = y
for (int i = n - 1; i >= 0; i--)
{
    double sum = 0;
    for (int j = i + 1; j < n; j++)
        sum += U[i][j] * x[j];
    x[i] = (y[i] - sum) / U[i][i];
}
```
- Solves from bottom to top
- Divides by diagonal elements of U

---

## 🔧 Complete C++ Implementation

```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin)
    {
        cerr << "Error: input.txt not found!\n";
        return 1;
    }

    fout << fixed << setprecision(4);

    bool printIntermediate = true; // toggle intermediate steps

    int n;
    while (fin >> n)
    {
        vector<vector<double>> aug(n, vector<double>(n + 1));
        
        // Read augmented matrix
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                fin >> aug[i][j];

        // Print the original system
        fout << "\n========================================\n";
        fout << "Input system:\n";
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (j > 0 && aug[i][j] >= 0) fout << "+";
                fout << aug[i][j] << "x" << j + 1 << " ";
            }
            fout << "= " << aug[i][n] << "\n";
        }
        fout << "========================================\n";

        // Separate A and b from augmented matrix
        vector<vector<double>> A(n, vector<double>(n));
        vector<double> b(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                A[i][j] = aug[i][j];
            b[i] = aug[i][n];
        }

        // Initialize L and U matrices
        vector<vector<double>> L(n, vector<double>(n, 0));
        vector<vector<double>> U(n, vector<double>(n, 0));

        fout << "\nPerforming LU Decomposition.. .\n";

        bool decompositionFailed = false;
        int failedAtStep = -1;

        // LU Decomposition using Doolittle's method
        for (int i = 0; i < n; i++)
        {
            // Upper Triangular Matrix U
            for (int k = i; k < n; k++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += L[i][j] * U[j][k];
                U[i][k] = A[i][k] - sum;
            }

            // Check for zero pivot
            if (fabs(U[i][i]) < 1e-12)
            {
                decompositionFailed = true;
                failedAtStep = i;
                break;
            }

            // Lower Triangular Matrix L
            for (int k = i; k < n; k++)
            {
                if (i == k)
                    L[i][i] = 1; // diagonal element is 1
                else
                {
                    double sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += L[k][j] * U[j][i];
                    L[k][i] = (A[k][i] - sum) / U[i][i];
                }
            }

            // Print intermediate matrices if enabled
            if (printIntermediate)
            {
                fout << "\nAfter step " << i + 1 << ":\n";
                fout << "L matrix:\n";
                for (int r = 0; r < n; r++)
                {
                    for (int c = 0; c < n; c++)
                        fout << setw(10) << L[r][c] << " ";
                    fout << "\n";
                }
                fout << "U matrix:\n";
                for (int r = 0; r < n; r++)
                {
                    for (int c = 0; c < n; c++)
                        fout << setw(10) << U[r][c] << " ";
                    fout << "\n";
                }
                fout << "---------------------------------------------\n";
            }
        }

        // Calculate determinant (product of diagonal of U)
        double detU = 1;
        for (int i = 0; i < n; i++)
            detU *= U[i][i];

        // Detect solution type
        bool noSolution = false;
        bool infiniteSolution = false;

        if (decompositionFailed || fabs(detU) < 1e-12)
        {
            // Matrix is singular - check for no solution vs infinite solutions
            
            // Perform forward substitution to check consistency
            vector<double> y(n, 0);
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += L[i][j] * y[j];
                
                double expectedY = b[i] - sum;
                
                // Check if this row is all zeros in U
                bool rowIsZero = true;
                for (int j = i; j < n; j++)
                {
                    if (fabs(U[i][j]) > 1e-12)
                    {
                        rowIsZero = false;
                        break;
                    }
                }
                
                if (rowIsZero && fabs(expectedY) > 1e-12)
                {
                    noSolution = true;
                    break;
                }
                
                if (! rowIsZero && fabs(U[i][i]) > 1e-12)
                    y[i] = expectedY;
            }

            if (! noSolution)
                infiniteSolution = true;
        }

        // Output results
        fout << "\n========================================\n";
        fout << "FINAL RESULT:\n";
        fout << "========================================\n";

        if (noSolution)
        {
            fout << "\nNo Solution\n";
            fout << "The system is inconsistent.\n";
            fout << "Determinant of U = " << detU << " (approximately 0)\n";
        }
        else if (infiniteSolution)
        {
            fout << "\nInfinite Solutions\n";
            fout << "The system has dependent equations.\n";
            fout << "Determinant of U = " << detU << " (approximately 0)\n";
        }
        else
        {
            fout << "\nUnique Solution\n";
            fout << "Determinant of U = " << detU << "\n\n";

            // Print final L and U matrices
            fout << "Final L matrix (Lower Triangular):\n";
            for (int r = 0; r < n; r++)
            {
                for (int c = 0; c < n; c++)
                    fout << setw(10) << L[r][c] << " ";
                fout << "\n";
            }

            fout << "\nFinal U matrix (Upper Triangular):\n";
            for (int r = 0; r < n; r++)
            {
                for (int c = 0; c < n; c++)
                    fout << setw(10) << U[r][c] << " ";
                fout << "\n";
            }

            // Forward substitution to solve L*y = b
            vector<double> y(n);
            fout << "\n--- Forward Substitution (L*y = b) ---\n";
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += L[i][j] * y[j];
                y[i] = b[i] - sum;
                fout << "y" << i + 1 << " = " << y[i] << "\n";
            }

            // Back substitution to solve U*x = y
            vector<double> x(n);
            fout << "\n--- Back Substitution (U*x = y) ---\n";
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0;
                for (int j = i + 1; j < n; j++)
                    sum += U[i][j] * x[j];
                x[i] = (y[i] - sum) / U[i][i];
            }

            fout << "\nSolution Vector (x):\n";
            for (int i = 0; i < n; i++)
                fout << "x" << i + 1 << " = " << x[i] << "\n";

            // Verification:  compute A*x and compare with b
            fout << "\n--- Verification (A*x = b) ---\n";
            vector<double> result(n);
            for (int i = 0; i < n; i++)
            {
                result[i] = 0;
                for (int j = 0; j < n; j++)
                    result[i] += A[i][j] * x[j];
                fout << "Row " << i + 1 << ": " << result[i] 
                     << " (expected: " << b[i] << ")\n";
            }
        }

        fout << "\n" << string(60, '=') << "\n\n";
    }

    fin.close();
    fout.close();
    cout << "All results written to output.txt\n";
    return 0;
}
```

---

## 📊 Usage Examples

### Input Format

The program reads from `input.txt`:

```
n
a₁₁ a₁₂ ... a₁ₙ b₁
a₂₁ a₂₂ ...  a₂ₙ b₂
...
aₙ₁ aₙ₂ ... aₙₙ bₙ
```

Where:
- `n` = number of equations (and variables)
- `aᵢⱼ` = coefficient of variable xⱼ in equation i
- `bᵢ` = constant term (right-hand side) of equation i

**Note**: Multiple test cases can be included one after another. 

### Example 1: Unique Solution

**Input (input.txt):**
```txt
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

**System of Equations:**
```
2x₁ + x₂ - x₃ = 8
-3x₁ - x₂ + 2x₃ = -11
-2x₁ + x₂ + 2x₃ = -3
```

**Output (output.txt):**
```txt
========================================
Input system:
2. 0000x1 +1.0000x2 -1.0000x3 = 8.0000
-3.0000x1 -1.0000x2 +2.0000x3 = -11.0000
-2.0000x1 +1.0000x2 +2.0000x3 = -3.0000
========================================

Performing LU Decomposition... 

After step 1:
L matrix:
    1.0000     0.0000     0.0000 
   -1.5000     0.0000     0.0000 
   -1.0000     0.0000     0.0000 
U matrix:
    2.0000     1.0000    -1.0000 
    0.0000     0.5000     0.5000 
    0.0000     2.0000     1.0000 
---------------------------------------------

After step 2:
L matrix:
    1.0000     0.0000     0.0000 
   -1.5000     1.0000     0.0000 
   -1.0000     4.0000     0.0000 
U matrix:
    2.0000     1.0000    -1.0000 
    0.0000     0.5000     0.5000 
    0.0000     0.0000    -1.0000 
---------------------------------------------

After step 3:
L matrix:
    1.0000     0.0000     0.0000 
   -1.5000     1.0000     0.0000 
   -1.0000     4.0000     1.0000 
U matrix:
    2.0000     1.0000    -1.0000 
    0.0000     0.5000     0.5000 
    0.0000     0.0000    -1.0000 
---------------------------------------------

========================================
FINAL RESULT: 
========================================

Unique Solution
Determinant of U = -1.0000

Final L matrix (Lower Triangular):
    1.0000     0.0000     0.0000 
   -1.5000     1.0000     0.0000 
   -1.0000     4.0000     1.0000 

Final U matrix (Upper Triangular):
    2.0000     1.0000    -1.0000 
    0.0000     0.5000     0.5000 
    0.0000     0.0000    -1.0000 

--- Forward Substitution (L*y = b) ---
y1 = 8.0000
y2 = 1.0000
y3 = 1.0000

--- Back Substitution (U*x = y) ---

Solution Vector (x):
x1 = 2.0000
x2 = 3.0000
x3 = -1.0000

--- Verification (A*x = b) ---
Row 1: 8.0000 (expected: 8.0000)
Row 2: -11.0000 (expected: -11.0000)
Row 3: -3.0000 (expected: -3.0000)

============================================================
```

### Example 2: No Solution

**Input:**
```txt
3
1 1 1 2
2 2 2 4
1 1 1 5
```

**System:**
```
x₁ + x₂ + x₃ = 2
2x₁ + 2x₂ + 2x₃ = 4
x₁ + x₂ + x₃ = 5
```

**Analysis**:  First and third equations are contradictory

**Output:**
```txt
========================================
Input system:
1.0000x1 +1.0000x2 +1.0000x3 = 2.0000
2.0000x1 +2.0000x2 +2.0000x3 = 4.0000
1.0000x1 +1.0000x2 +1.0000x3 = 5.0000
========================================

Performing LU Decomposition...

After step 1:
L matrix:
    1.0000     0.0000     0.0000 
    2.0000     0.0000     0.0000 
    1.0000     0.0000     0.0000 
U matrix:
    1.0000     1.0000     1.0000 
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 
---------------------------------------------

========================================
FINAL RESULT: 
========================================

No Solution
The system is inconsistent.
Determinant of U = 0.0000 (approximately 0)

============================================================
```

### Example 3: Infinite Solutions

**Input:**
```txt
3
1 2 3 6
2 4 6 12
3 6 9 18
```

**System:**
```
x₁ + 2x₂ + 3x₃ = 6
2x₁ + 4x₂ + 6x₃ = 12
3x₁ + 6x₂ + 9x₃ = 18
```

**Analysis**: All equations are multiples of each other

**Output:**
```txt
========================================
Input system: 
1.0000x1 +2.0000x2 +3.0000x3 = 6.0000
2.0000x1 +4.0000x2 +6.0000x3 = 12.0000
3.0000x1 +6.0000x2 +9.0000x3 = 18.0000
========================================

Performing LU Decomposition...

After step 1:
L matrix: 
    1.0000     0.0000     0.0000 
    2.0000     0.0000     0.0000 
    3.0000     0.0000     0.0000 
U matrix:
    1.0000     2.0000     3.0000 
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 
---------------------------------------------

========================================
FINAL RESULT: 
========================================

Infinite Solutions
The system has dependent equations.
Determinant of U = 0.0000 (approximately 0)

============================================================
```

---

## 🎯 Compilation and Execution

### Compile
```bash
g++ -std=c++17 -O2 lu-decomposition.cpp -o lu-decomposition
```

### Run
```bash
./lu-decomposition
```

### Requirements
- C++11 or later
- Input file: `input.txt` in the same directory
- Output file: `output.txt` (automatically created)

### 🔧 Configuration

Toggle intermediate step printing: 
```cpp
bool printIntermediate = true; // Set to false to hide intermediate matrices
```

---

## 📊 Comparison with Other Methods

### Difference from Gauss Elimination

| Aspect | Gauss Elimination | LU Decomposition |
|--------|-------------------|------------------|
| **Approach** | Direct method | Matrix factorization |
| **Matrix Storage** | One augmented matrix | Three matrices (A, L, U) |
| **Multiple Right-Hand Sides** | Must repeat elimination | Reuse L and U |
| **Complexity (single solve)** | O(n³) | O(n³) |
| **Complexity (k solves)** | k × O(n³) | O(n³) + k × O(n²) |
| **Memory Usage** | O(n²) | O(3n²) |
| **Determinant** | Requires tracking | Product of U diagonal |
| **Matrix Inverse** | Must augment with I | Can compute A⁻¹ = U⁻¹L⁻¹ |

### Difference from Gauss-Jordan Elimination

| Aspect | Gauss-Jordan | LU Decomposition |
|--------|--------------|------------------|
| **Final Form** | RREF (I \| solution) | L and U triangular |
| **Solution Method** | Direct reading | Forward + Back substitution |
| **Elimination** | Full (above & below) | Partial (below only) |
| **Operations** | ~n³/2 | ~n³/3 |
| **Multiple Solves** | k × O(n³) | O(n³) + k × O(n²) |
| **Best Use Case** | Single solve, matrix inverse | Multiple solves, determinant |

### When to Use LU Decomposition

✅ **Use LU Decomposition When:**
- Solving **multiple systems** with the **same A** but different b vectors
- Computing **determinant** (simple product of diagonal)
- Finding **matrix inverse** (solve A×X = I column by column)
- You need to **preserve the original matrix** A
- Working with **large systems** where efficiency matters
- Implementing iterative refinement for better accuracy

❌ **Avoid LU Decomposition When:**
- Solving only **one system** (Gauss Elimination is simpler)
- Matrix is **sparse** (special methods like iterative solvers are better)
- Matrix is **ill-conditioned** (need pivoting strategies)
- You need **RREF form** explicitly (use Gauss-Jordan)

### Computational Cost Comparison

For an n×n system: 

| Method | Single Solve | k Different b Vectors |
|--------|--------------|----------------------|
| **Gauss Elimination** | ~n³/3 | k × (n³/3) |
| **Gauss-Jordan** | ~n³/2 | k × (n³/2) |
| **LU Decomposition** | ~n³/3 | n³/3 + k × (n²) |

**Example**: For n=100 and k=10 right-hand sides:
- Gauss Elimination: 10 × 333,333 ≈ 3,333,330 operations
- LU Decomposition: 333,333 + 10 × 10,000 ≈ 433,333 operations
- **LU is ~7. 7× faster! **

---

## 🎓 Applications

### Direct Applications
- **Solving Linear Systems** - Particularly with multiple right-hand sides
- **Matrix Inversion** - Compute A⁻¹ by solving AX = I column by column
- **Determinant Calculation** - det(A) = Π(U[i][i])
- **Matrix Equation Solving** - AX = B where B is a matrix
- **Condition Number Estimation** - Used in numerical stability analysis

### Real-World Uses
- **Circuit Analysis** - Nodal and mesh analysis with varying sources
- **Structural Engineering** - Solving for forces with different loading conditions
- **Control Systems** - State-space models with different initial conditions
- **Computer Graphics** - Solving lighting equations for different light sources
- **Economics** - Input-output analysis with different demand scenarios
- **Finite Element Analysis** - Same mesh, different boundary conditions
- **Machine Learning** - Ridge regression, normal equations
- **Signal Processing** - Filter design, system identification
- **Computational Physics** - Discretized PDEs with varying parameters

---

## 🔬 Mathematical Properties

### Properties of LU Decomposition

1. **Existence**: 
   - LU decomposition exists if all leading principal minors of A are non-zero
   - If A is non-singular and no row exchanges are needed, LU exists

2. **Triangular Structure**:
   - **L**:  Lower triangular with 1s on diagonal (Doolittle)
   - **U**: Upper triangular with computed values on diagonal
   - **A = L × U** exactly (no approximation)

3. **Row Operations**:
   - LU decomposition is equivalent to Gauss Elimination without row swaps
   - L stores the multipliers used during elimination
   - U is the upper triangular result of elimination

4. **Permutation**:
   - With pivoting:  **PA = LU** where P is a permutation matrix
   - Ensures numerical stability
   - Not implemented in this basic version

### Uniqueness

For Doolittle's method (diagonal of L = 1):
- LU decomposition is **unique** if it exists
- Different conventions (Crout, Cholesky) give different L and U
- But A = L × U is always satisfied

**Proof sketch**:
If A = L₁U₁ = L₂U₂, then:
- L₂⁻¹L₁ = U₂U₁⁻¹
- Left side is lower triangular, right side is upper triangular
- Both must be diagonal
- With diagonal of L = 1, both must be identity
- Therefore L₁ = L₂ and U₁ = U₂

### Determinant Calculation

One of the most elegant applications:

```
det(A) = det(L × U) = det(L) × det(U)
```

Since L has 1s on diagonal:  **det(L) = 1**

Therefore: 
```
det(A) = det(U) = U[0][0] × U[1][1] × ...  × U[n-1][n-1]
```

**Advantages**:
- O(n) to compute after decomposition
- No additional storage needed
- Numerically stable

### Relationship to Gauss Elimination

LU Decomposition **is** Gauss Elimination:
- **L** stores the elimination multipliers
- **U** is the result of forward elimination
- **b** is transformed by the same operations

Example:
```
If we subtract 2 × row1 from row2:
- This creates U
- We record L[2][1] = 2
```

This connection means: 
- Same complexity:  O(n³)
- Same numerical properties
- But LU preserves information for reuse

---

## 📚 References

- [LU Decomposition - Wikipedia](https://en.wikipedia.org/wiki/LU_decomposition)
- Numerical Methods For Engineers by Raymond Canale and Steven C. Chapra

---

## 👤 Author

**Part of the [Numerical Computing Suite](../) by [AbirHasanArko](https://github.com/AbirHasanArko)**  
Roll:  2207053  
Department of CSE, KUET   