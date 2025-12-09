# Gauss-Jordan Elimination Method

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](gauss-jordan-elimination.cpp)
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
  - [Solution Detection Logic](#solution-detection-logic)
- [Complete C++ Implementation](#-complete-c-implementation)
- [Usage Examples](#-usage-examples)
  - [Example 1: Unique Solution](#example-1-unique-solution)
  - [Example 2: No Solution](#example-2-no-solution)
  - [Example 3: Infinite Solutions](#example-3-infinite-solutions)
- [Compilation and Execution](#-compilation-and-execution)
- [Comparison:  When to Use Which Method](#-comparison-when-to-use-which-method)
  - [Use Gauss-Jordan When](#use-gauss-jordan-when)
  - [Use Gauss Elimination When](#use-gauss-elimination-when)
  - [Computational Cost Comparison](#computational-cost-comparison)
- [Applications](#-applications)
  - [Direct Applications](#direct-applications)
  - [Real-World Uses](#real-world-uses)
- [Mathematical Properties](#-mathematical-properties)
  - [Properties of RREF](#properties-of-rref)
  - [Uniqueness](#uniqueness)
  - [Relationship to Matrix Inverse](#relationship-to-matrix-inverse)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

The **Gauss-Jordan Elimination Method** is an extension of Gaussian Elimination that transforms the coefficient matrix into **reduced row echelon form (RREF)** rather than just upper triangular form. This method eliminates the need for back substitution, as the solution can be read directly from the final matrix.

This implementation in C++ solves systems of linear equations and can detect three types of solutions:
- **Unique Solution** - The system has exactly one solution
- **No Solution** - The system is inconsistent
- **Infinite Solutions** - The system has infinitely many solutions

### ⚙️ Features

✅ **Partial Pivoting** - Selects best pivot for numerical stability  
✅ **Row Normalization** - Makes diagonal elements equal to 1  
✅ **Full Elimination** - Creates zeros above and below pivots  
✅ **Solution Detection** - Identifies unique, no, or infinite solutions  
✅ **Direct Reading** - No back substitution required  
✅ **Multiple Test Cases** - Process several systems in one run  
✅ **Intermediate Output** - View matrix transformation steps  
✅ **File I/O** - Reads from input.txt, writes to output.txt  
✅ **High Precision** - Uses double precision with 2 decimal places  
✅ **Error Handling** - Checks for missing input file  

---

## 🧮 Theory & Algorithm

### Mathematical Foundation

Given a system of linear equations:
```
a₁₁x₁ + a₁₂x₂ + ... + a₁ₙxₙ = b₁
a₂₁x₁ + a₂₂x₂ + ... + a₂ₙxₙ = b₂
...  
aₙ₁x₁ + aₙ₂x₂ + ...  + aₙₙxₙ = bₙ
```

The method converts the augmented matrix [A|b] into reduced row echelon form where the coefficient matrix becomes an identity matrix. 

### Algorithm Steps

1. **Forward Elimination Phase:**
   - For each pivot position (i, i):
     - **Partial Pivoting**: Find the row with the maximum absolute value in column i (from row i to n)
     - Swap the current row with the maximum row
     - **Normalize the pivot row**:  Divide the entire row by the pivot element to make pivot = 1
     - Eliminate all entries above AND below the pivot by subtracting multiples of the pivot row
     - Calculate the multiplication factor:  `factor = a[k][i]` (since pivot is now 1)
     - Update row k: `a[k][j] = a[k][j] - factor × a[i][j]` for all rows k ≠ i

2. **Solution Detection Phase:**
   - Calculate the rank of the coefficient matrix
   - Check for inconsistency:  If any row has all zeros in coefficients but non-zero constant → **No Solution**
   - Check for infinite solutions: If rank < n → **Infinite Solutions**
   - Otherwise → **Unique Solution**, read directly from the matrix

3. **Direct Solution Reading** (for unique solutions):
   - No back substitution needed! 
   - The matrix is already in reduced row echelon form (RREF)
   - Simply read:  `xᵢ = a[i][n]` for i = 1 to n

### Complexity Analysis

#### Time Complexity
- **Forward Elimination**: O(n³)
- **No Back Substitution**: O(n) - just reading values
- **Overall**: O(n³)

#### Space Complexity
- O(n²) for storing the augmented matrix

---

## 💻 Implementation Details

### Key Components

```cpp
vector<vector<double>> a(n, vector<double>(n + 1));
```
- Creates an augmented matrix of size n×(n+1) to store coefficients and constants

```cpp
// Partial Pivoting
int maxRow = i;
for (int k = i + 1; k < n; k++)
    if (fabs(a[k][i]) > fabs(a[maxRow][i]))
        maxRow = k;
swap(a[i], a[maxRow]);
```
- Finds the row with maximum absolute value in current column
- Swaps it to pivot position for numerical stability

```cpp
// Make diagonal element 1 (Row Normalization)
double pivot = a[i][i];
for (int j = i; j <= n; j++)
    a[i][j] /= pivot;
```
- **Key Step**: Divides the entire pivot row by the pivot element
- Ensures diagonal element becomes exactly 1
- This is the main difference from standard Gauss Elimination

```cpp
// Eliminate entries ABOVE AND BELOW the pivot
for (int k = 0; k < n; k++)
{
    if (k != i)
    {
        double factor = a[k][i];
        for (int j = i; j <= n; j++)
            a[k][j] -= factor * a[i][j];
    }
}
```
- **Critical Difference**: Eliminates ALL rows (not just below)
- Makes all entries above and below the pivot equal to zero
- Creates the reduced row echelon form

```cpp
for (int i = 0; i < n; i++)
    fout << "x" << i + 1 << " = " << a[i][n] << "\n";
```
- Direct solution reading - no back substitution loop needed
- The solution is already available in the last column of the RREF matrix

#### Solution Detection Logic

```cpp
bool noSolution = false, infiniteSolution = false;
int rank = 0;

for (int i = 0; i < n; i++)
{
    bool allZero = true;
    for (int j = 0; j < n; j++)
        if (fabs(a[i][j]) > 1e-12) allZero = false;

    if (allZero && fabs(a[i][n]) > 1e-12)
    {
        noSolution = true;
        break;
    }
    if (!allZero)
        rank++;
}
```
- Checks each row for consistency
- A row like `[0 0 0 | 5]` indicates no solution (0 = 5 is impossible)
- Counts non-zero rows to determine rank

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

    fout << fixed << setprecision(2);

    bool printIntermediate = true; // toggle intermediate steps

    int n;
    while (fin >> n)
    {
        vector<vector<double>> a(n, vector<double>(n + 1));
        vector<double> x(n);

        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                fin >> a[i][j];

        // Print the original system
        fout << "\nInput system:\n";
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (j > 0 && a[i][j] >= 0) fout << "+";
                fout << a[i][j] << "x" << j + 1 << " ";
            }
            fout << "= " << a[i][n] << "\n";
        }

        // Gauss-Jordan Elimination
        int step = 1;
        for (int i = 0; i < n; i++)
        {
            // Partial Pivoting
            int maxRow = i;
            for (int k = i + 1; k < n; k++)
                if (fabs(a[k][i]) > fabs(a[maxRow][i]))
                    maxRow = k;
            
            if (i != maxRow)
            {
                swap(a[i], a[maxRow]);
                if (printIntermediate)
                {
                    fout << "\nAfter swapping row " << i + 1 << " with row " << maxRow + 1 << ":\n";
                    for (int r = 0; r < n; r++)
                    {
                        for (int c = 0; c <= n; c++)
                            fout << a[r][c] << "\t";
                        fout << "\n";
                    }
                }
            }

            if (fabs(a[i][i]) < 1e-12)
                continue;

            // Make diagonal element 1 (normalize pivot row)
            double pivot = a[i][i];
            for (int j = i; j <= n; j++)
                a[i][j] /= pivot;

            if (printIntermediate)
            {
                fout << "\nStep " << step++ << " - Making diagonal element a[" << i + 1 << "][" << i + 1 << "] = 1:\n";
                for (int r = 0; r < n; r++)
                {
                    for (int c = 0; c <= n; c++)
                        fout << a[r][c] << "\t";
                    fout << "\n";
                }
            }

            // Eliminate column i in ALL other rows (both above and below)
            for (int k = 0; k < n; k++)
            {
                if (k != i)
                {
                    double factor = a[k][i];
                    for (int j = i; j <= n; j++)
                        a[k][j] -= factor * a[i][j];
                }
            }

            // Print intermediate matrix if enabled
            if (printIntermediate)
            {
                fout << "\nStep " << step++ << " - Eliminating column " << i + 1 << " in all other rows:\n";
                for (int r = 0; r < n; r++)
                {
                    for (int c = 0; c <= n; c++)
                        fout << a[r][c] << "\t";
                    fout << "\n";
                }
            }
        }

        // Detect solution type
        bool noSolution = false, infiniteSolution = false;
        int rank = 0;

        for (int i = 0; i < n; i++)
        {
            bool allZero = true;
            for (int j = 0; j < n; j++)
                if (fabs(a[i][j]) > 1e-12) allZero = false;

            if (allZero && fabs(a[i][n]) > 1e-12)
            {
                noSolution = true;
                break;
            }
            if (! allZero)
                rank++;
        }

        if (! noSolution && rank < n)
            infiniteSolution = true;

        // Write output
        fout << "\n========================================\n";
        fout << "FINAL RESULT:\n";
        fout << "========================================\n";

        if (noSolution)
        {
            fout << "\nNo Solution\n";
            fout << "The system is inconsistent.\n";
        }
        else if (infiniteSolution)
        {
            fout << "\nInfinite Solutions\n";
            fout << "The system has dependent equations.\n";
        }
        else
        {
            fout << "\nUnique Solution\n";
            fout << "\nFinal Reduced Row Echelon Form (RREF):\n";
            for (int r = 0; r < n; r++)
            {
                for (int c = 0; c <= n; c++)
                    fout << a[r][c] << "\t";
                fout << "\n";
            }

            fout << "\nSolution:\n";
            for (int i = 0; i < n; i++)
            {
                x[i] = a[i][n];
                fout << "x" << i + 1 << " = " << x[i] << "\n";
            }
        }

        fout << "\n" << string(60, '=') << "\n\n";
    }

    fin.close();
    fout. close();
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

## 📤 Output Format

The program writes to `output.txt`:

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
Input system:
2.00x1 +1.00x2 -1.00x3 = 8.00
-3.00x1 -1.00x2 +2.00x3 = -11.00
-2.00x1 +1.00x2 +2.00x3 = -3.00

After swapping row 1 with row 2:
-3.00	-1.00	2.00	-11.00	
2.00	1.00	-1.00	8.00	
-2.00	1.00	2.00	-3.00	

Step 1 - Making diagonal element a[1][1] = 1:
1.00	0.33	-0.67	3.67	
2.00	1.00	-1.00	8.00	
-2.00	1.00	2.00	-3.00	

Step 2 - Eliminating column 1 in all other rows:
1.00	0.33	-0.67	3.67	
0.00	0.33	0.33	0.67	
0.00	1.67	0.67	4.33	

After swapping row 2 with row 3:
1.00	0.33	-0.67	3.67	
0.00	1.67	0.67	4.33	
0.00	0.33	0.33	0.67	

Step 3 - Making diagonal element a[2][2] = 1:
1.00	0.33	-0.67	3.67	
0.00	1.00	0.40	2.60	
0.00	0.33	0.33	0.67	

Step 4 - Eliminating column 2 in all other rows:
1.00	0.00	-0.80	2.80	
0.00	1.00	0.40	2.60	
0.00	0.00	0.20	-0.20	

Step 5 - Making diagonal element a[3][3] = 1:
1.00	0.00	-0.80	2.80	
0.00	1.00	0.40	2.60	
0.00	0.00	1.00	-1.00	

Step 6 - Eliminating column 3 in all other rows:
1.00	0.00	0.00	2.00	
0.00	1.00	0.00	3.00	
0.00	0.00	1.00	-1.00	

========================================
FINAL RESULT:
========================================

Unique Solution

Final Reduced Row Echelon Form (RREF):
1.00	0.00	0.00	2.00	
0.00	1.00	0.00	3.00	
0.00	0.00	1.00	-1.00	

Solution:
x1 = 2.00
x2 = 3.00
x3 = -1.00
```

### Example 2: No Solution

**Input:**
```
3
1 1 1 2
2 2 2 4
1 1 1 5
```

**System:**
```
x₁ + x₂ + x₃ = 2
2x₁ + 2x₂ + 2x₃ = 5
x₁ + x₂ + x₃ = 5
```

**Analysis**:  First and third equations are contradictory (same left side, different right side)

**Output:**
```
Input system:
1.00x1 +1.00x2 +1.00x3 = 2.00
2.00x1 +2.00x2 +2.00x3 = 4.00
1.00x1 +1.00x2 +1.00x3 = 5.00

After swapping row 1 with row 2:
2.00	2.00	2.00	4.00	
1.00	1.00	1.00	2.00	
1.00	1.00	1.00	5.00	

Step 1 - Making diagonal element a[1][1] = 1:
1.00	1.00	1.00	2.00	
1.00	1.00	1.00	2.00	
1.00	1.00	1.00	5.00	

Step 2 - Eliminating column 1 in all other rows:
1.00	1.00	1.00	2.00	
0.00	0.00	0.00	0.00	
0.00	0.00	0.00	3.00	

========================================
FINAL RESULT:
========================================

No Solution
The system is inconsistent.
```

### Example 3: Infinite Solutions

**Input:**
```
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

**Analysis**: All equations are multiples of each other (rank = 1 < 3)

**Output:**
```
Input system:
1.00x1 +2.00x2 +3.00x3 = 6.00
2.00x1 +4.00x2 +6.00x3 = 12.00
3.00x1 +6.00x2 +9.00x3 = 18.00

After swapping row 1 with row 3:
3.00	6.00	9.00	18.00	
2.00	4.00	6.00	12.00	
1.00	2.00	3.00	6.00	

Step 1 - Making diagonal element a[1][1] = 1:
1.00	2.00	3.00	6.00	
2.00	4.00	6.00	12.00	
1.00	2.00	3.00	6.00	

Step 2 - Eliminating column 1 in all other rows:
1.00	2.00	3.00	6.00	
0.00	0.00	0.00	0.00	
0.00	0.00	0.00	0.00	

========================================
FINAL RESULT:
========================================

Infinite Solutions
The system has dependent equations.
```

---

## 🎯 Compilation and Execution

### Compile
```bash
g++ -std=c++17 -O2 gauss-jordan-elimination-method.cpp -o gauss-jordan
```

### Run
```bash
./gauss-jordan
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

## 📊 Comparison: When to Use Which Method

### Difference from Gauss Elimination

| Aspect | Gauss Elimination | Gauss-Jordan Elimination |
|--------|-------------------|--------------------------|
| **Final Matrix Form** | Upper Triangular | Reduced Row Echelon Form (RREF) |
| **Diagonal Elements** | Non-zero values | All equal to 1 |
| **Off-diagonal Elements** | Upper triangle:  any value<br>Lower triangle: zero | All zeros |
| **Elimination Direction** | Below pivot only | Above AND below pivot |
| **Row Normalization** | Not required | Required (divide by pivot) |
| **Solution Method** | Back substitution | Direct reading |
| **Number of Operations** | ~n³/3 multiplications | ~n³/2 multiplications |
| **When to Use** | Just solving equations | Solving equations + Matrix inversion |

### Use Gauss-Jordan When: 
- We need the **inverse of a matrix** (augment with identity matrix)
- We want to **understand RREF** for educational purposes
- We need to **find the rank** explicitly
- We prefer a **simpler final step** (no back substitution)
- The system is **small to medium** sized (n < 1000)

### Use Gauss Elimination When:
- We only need to **solve the system** (not find inverse)
- **Performance is critical** (large systems)
- We want **fewer arithmetic operations**
- The system is **very large** (n > 1000)

### Computational Cost Comparison
For an n×n system: 

| Method | Multiplications/Divisions | Additions/Subtractions |
|--------|---------------------------|------------------------|
| **Gauss Elimination** | ~n³/3 + n²/2 | ~n³/3 + n²/2 - n |
| **Gauss-Jordan** | ~n³/2 + n² | ~n³/2 + n² - n |

Gauss-Jordan performs approximately **50% more operations**, but asymptotic complexity remains O(n³).

---

## 🎓 Applications

### Direct Applications
- **Solving Linear Systems** - Find unique solutions to Ax = b
- **Matrix Inversion** - Compute A⁻¹ by augmenting [A|I]
- **Rank Determination** - Number of non-zero rows in RREF
- **Testing Linear Independence** - Check if columns are independent

### Real-World Uses
- **Circuit Analysis** - Solving Kirchhoff's laws
- **Computer Graphics** - Transformation matrices
- **Economics** - Input-output models (Leontief models)
- **Engineering** - Structural analysis, heat transfer
- **Machine Learning** - Linear regression, least squares
- **Operations Research** - Linear programming (simplex method)
- **Chemistry** - Balancing chemical equations
- **Physics** - Systems of forces, electrical networks

---

## 🔬 Mathematical Properties

### Properties of RREF
1. **Leading Entry**: Each non-zero row has a leading 1
2. **Column Zeros**: All entries in a pivot column (except the pivot) are 0
3. **Staircase Pattern**: Each leading 1 is to the right of the leading 1 in the row above
4. **Zero Rows**: All zero rows are at the bottom

### Uniqueness
- The RREF of a matrix is **unique**
- Two different row operations sequences lead to the same RREF
- This property makes RREF useful for comparing matrices

### Relationship to Matrix Inverse
For an n×n matrix A:
- Augment:  [A | I]
- Apply Gauss-Jordan: [I | A⁻¹]
- If you can't reach [I | ? ], then A is singular (no inverse)

---

## 📚 References

- [Gauss-Jordan Elimination Method - GeeksforGeeks](https://www.geeksforgeeks.org/dsa/program-for-gauss-jordan-elimination-method/)
- Numerical Methods For Engineers by Raymond Canale and Steven C. Chapra

---

## 👤 Author

**Part of the [Numerical Computing Suite](../) by [AbirHasanArko](https://github.com/AbirHasanArko)**  
Roll: 2207053   
Department of CSE, KUET   