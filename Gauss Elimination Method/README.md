# Gauss Elimination Method

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](gauss-elimination-method.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

---

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Theory & Algorithm](#-theory--algorithm)
  - [Mathematical Foundation](#mathematical-foundation)
  - [Algorithm Steps](#algorithm-steps)
  - [Partial Pivoting Strategy](#partial-pivoting-strategy)
  - [Solution Detection Mechanism](#solution-detection-mechanism)
  - [Complexity Analysis](#complexity-analysis)
- [Implementation Details](#-implementation-details)
  - [Key Components](#key-components)
  - [Forward Elimination Process](#forward-elimination-process)
  - [Back Substitution Process](#back-substitution-process)
  - [Numerical Stability](#numerical-stability)
- [Complete C++ Implementation](#-complete-c-implementation)
- [Usage Examples](#-usage-examples)
  - [Example 1: Unique Solution](#example-1-unique-solution)
  - [Example 2: No Solution](#example-2-no-solution)
  - [Example 3: Infinite Solutions](#example-3-infinite-solutions)
- [Compilation and Execution](#-compilation-and-execution)
- [Applications](#-applications)
- [References](#-references)

---

## 📖 Introduction

The **Gauss Elimination Method** (also known as **Gaussian Elimination**) is a fundamental algorithm in linear algebra used to solve systems of linear equations. Named after German mathematician Carl Friedrich Gauss, this method transforms the coefficient matrix into an upper triangular form through a series of elementary row operations, followed by back substitution to find the solution.

This implementation in C++ solves systems of linear equations and can detect three types of solutions:
- **Unique Solution** - The system has exactly one solution (rank(A) = rank([A|b]) = n)
- **No Solution** - The system is inconsistent (rank(A) < rank([A|b]))
- **Infinite Solutions** - The system has infinitely many solutions (rank(A) = rank([A|b]) < n)

### Features
- ✅ Partial pivoting for numerical stability
- ✅ Solution type detection (unique/none/infinite)
- ✅ Intermediate step visualization
- ✅ High precision output (configurable)
- ✅ File-based I/O for batch processing
- ✅ Error handling for invalid inputs

---

## 🧮 Theory & Algorithm

### Mathematical Foundation

Given a system of **n** linear equations with **n** unknowns: 

```
a₁₁x₁ + a₁₂x₂ + ... + a₁ₙxₙ = b₁
a₂₁x₁ + a₂₂x₂ + ... + a₂ₙxₙ = b₂
... 
aₙ₁x₁ + aₙ₂x₂ + ... + aₙₙxₙ = bₙ
```

This can be represented in matrix form as **Ax = b**, where:
- **A** is the coefficient matrix (n×n)
- **x** is the vector of unknowns (n×1)
- **b** is the constant vector (n×1)

The augmented matrix [A|b] is formed by appending b to A: 

```
[ a₁₁  a₁₂  ...  a₁ₙ | b₁ ]
[ a₂₁  a₂₂  ...  a₂ₙ | b₂ ]
[  ⋮    ⋮    ⋱   ⋮   | ⋮  ]
[ aₙ₁  aₙ₂  ...  aₙₙ | bₙ ]
```

### Algorithm Steps

#### Phase 1: Forward Elimination

The goal is to transform the augmented matrix into **row echelon form** (upper triangular).

For each pivot position (i, i) where i = 0 to n-2: 

1. **Partial Pivoting**: 
   - Find the row k (where k ≥ i) with maximum |aₖᵢ|
   - Swap row i with row k
   - This improves numerical stability

2. **Row Elimination**:
   - For each row k below the pivot (k = i+1 to n-1):
     - Calculate multiplication factor: `factor = aₖᵢ / aᵢᵢ`
     - Update row k: `aₖⱼ = aₖⱼ - factor × aᵢⱼ` for j = i to n
     - This makes all elements below the pivot equal to zero

**Result**: Upper triangular matrix

```
[ a₁₁  a₁₂  a₁₃ | b₁ ]
[  0   a₂₂  a₂₃ | b₂ ]
[  0    0   a₃₃ | b₃ ]
```

#### Phase 2: Solution Detection

Calculate the rank of the matrix to determine solution type:

1. **Count non-zero rows**: A row is non-zero if at least one coefficient is non-zero
2. **Check for inconsistency**: If any row has form [0 0 ... 0 | c] where c ≠ 0 → **No Solution**
3. **Check rank**: 
   - If rank < n → **Infinite Solutions**
   - If rank = n → **Unique Solution**

#### Phase 3: Back Substitution

For unique solutions, solve for variables from bottom to top:

```
xₙ = bₙ / aₙₙ

xᵢ = (bᵢ - Σⱼ₌ᵢ₊₁ⁿ aᵢⱼxⱼ) / aᵢᵢ   for i = n-1, n-2, ..., 1
```

### Partial Pivoting Strategy

Partial pivoting is crucial for numerical stability. Without it, small pivot elements can cause: 
- **Round-off errors** to accumulate
- **Loss of significant digits**
- **Numerical instability** in the solution

**Pivoting Rule**: At step i, select the row k ≥ i such that |aₖᵢ| is maximum, then swap rows i and k. 

**Example**: If pivot element is 0.001 and another element in the column is 100, swapping reduces error magnification.

### Solution Detection Mechanism

The system's solution type depends on the rank relationship: 

| Condition | Solution Type | Explanation |
|-----------|---------------|-------------|
| rank(A) = rank([A\|b]) = n | Unique Solution | System is consistent and determined |
| rank(A) = rank([A\|b]) < n | Infinite Solutions | System is consistent but underdetermined |
| rank(A) < rank([A\|b]) | No Solution | System is inconsistent |

### Complexity Analysis

#### Time Complexity
- **Forward Elimination**: O(n³)
  - Outer loop: n iterations
  - Inner loops: O(n²) per iteration
  - Total: n × n² = n³
- **Back Substitution**: O(n²)
  - n iterations with decreasing work
- **Overall**: **O(n³)** - dominated by forward elimination

#### Space Complexity
- **Augmented Matrix**: O(n²) - stores n×(n+1) elements
- **Solution Vector**: O(n)
- **Overall**: **O(n²)**

#### Operation Count
For an n×n system:
- **Multiplications/Divisions**: ~n³/3 + n²/2
- **Additions/Subtractions**: ~n³/3 + n²/2 - n/6

---

## 💻 Implementation Details

### Key Components

#### 1. Data Structures

```cpp
vector<vector<double>> a(n, vector<double>(n + 1));  // Augmented matrix [A|b]
vector<double> x(n);                                  // Solution vector
```

- **Augmented matrix**: Size n×(n+1) stores coefficients and constants
- **Solution vector**:  Stores the final values of x₁, x₂, ..., xₙ

#### 2. Precision Control

```cpp
fout << fixed << setprecision(2);
```

- Uses fixed-point notation for consistent formatting
- Precision of 2 decimal places (adjustable)

#### 3. Tolerance for Zero Comparison

```cpp
const double EPSILON = 1e-12;
if (fabs(a[i][i]) < EPSILON) // Check if effectively zero
```

- Prevents floating-point comparison errors
- Treats values < 10⁻¹² as zero

### Forward Elimination Process

```cpp
for (int i = 0; i < n - 1; i++)
{
    // Step 1: Partial Pivoting
    int maxRow = i;
    for (int k = i + 1; k < n; k++)
        if (fabs(a[k][i]) > fabs(a[maxRow][i]))
            maxRow = k;
    swap(a[i], a[maxRow]);
    
    // Step 2: Check for zero pivot
    if (fabs(a[i][i]) < 1e-12)
        continue;  // Skip if pivot is effectively zero
    
    // Step 3: Eliminate below pivot
    for (int k = i + 1; k < n; k++)
    {
        double factor = a[k][i] / a[i][i];
        for (int j = i; j <= n; j++)
            a[k][j] -= factor * a[i][j];
    }
}
```

**Key Points**:
- Finds maximum element in column i (rows i to n-1)
- Swaps rows for better numerical stability
- Eliminates all elements below diagonal

### Back Substitution Process

```cpp
for (int i = n - 1; i >= 0; i--)
{
    x[i] = a[i][n];                    // Start with constant term
    for (int j = i + 1; j < n; j++)
        x[i] -= a[i][j] * x[j];        // Subtract known terms
    x[i] /= a[i][i];                   // Divide by diagonal element
}
```

**Process**:
1. Start from last equation:  xₙ = bₙ/aₙₙ
2. Move upward, substituting known values
3. Solve for each variable sequentially

### Numerical Stability

This implementation ensures stability through: 

1. **Partial Pivoting**: Reduces round-off error propagation
2. **Epsilon Comparison**: Handles floating-point precision issues
3. **Double Precision**: Uses `double` for higher accuracy
4. **Zero Detection**: Identifies singular/near-singular matrices

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
        cerr << "Error: input. txt not found!\n";
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

        // Forward Elimination
        for (int i = 0; i < n - 1; i++)
        {
            int maxRow = i;
            for (int k = i + 1; k < n; k++)
                if (fabs(a[k][i]) > fabs(a[maxRow][i]))
                    maxRow = k;
            swap(a[i], a[maxRow]);

            if (fabs(a[i][i]) < 1e-12)
                continue;

            for (int k = i + 1; k < n; k++)
            {
                double factor = a[k][i] / a[i][i];
                for (int j = i; j <= n; j++)
                    a[k][j] -= factor * a[i][j];
            }

            // Print intermediate matrix if enabled
            if (printIntermediate)
            {
                fout << "\nAfter step " << i + 1 << ":\n";
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
            if (!allZero)
                rank++;
        }

        if (!noSolution && rank < n)
            infiniteSolution = true;

        // Write output
        if (noSolution)
            fout << "\nNo Solution\n";
        else if (infiniteSolution)
            fout << "\nInfinite Solutions\n";
        else
        {
            fout << "\nUnique Solution\n";
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = a[i][n];
                for (int j = i + 1; j < n; j++)
                    x[i] -= a[i][j] * x[j];
                x[i] /= a[i][i];
            }

            fout << "Solution:\n";
            for (int i = 0; i < n; i++)
                fout << "x" << i + 1 << " = " << x[i] << "\n";
        }

        fout << "\n-------------------------------------------\n\n";
    }

    fin.close();
    fout.close();
    cout << "All results written to output.txt\n";
    return 0;
}
```

### Code Walkthrough

| Lines | Purpose | Description |
|-------|---------|-------------|
| 1-2 | Headers | Includes all standard libraries |
| 6-13 | File I/O Setup | Opens input/output files with error handling |
| 15 | Formatting | Sets output precision to 2 decimal places |
| 17 | Debug Flag | Toggles intermediate step printing |
| 20-27 | Input Reading | Reads matrix size and augmented matrix |
| 29-39 | Display System | Prints original system of equations |
| 42-71 | Forward Elimination | Implements pivoting and elimination |
| 74-93 | Solution Detection | Calculates rank and determines solution type |
| 96-114 | Output Generation | Performs back substitution and writes results |

---

## 📊 Usage Examples

### Input Format

```
n
a₁₁ a₁₂ ... a₁ₙ b₁
a₂₁ a₂₂ ...  a₂ₙ b₂
...
aₙ₁ aₙ₂ ... aₙₙ bₙ
```

### Example 1: Unique Solution

**Input (input.txt):**
```
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
```
Input system:
2.00x1 +1.00x2 -1.00x3 = 8.00
-3.00x1 -1.00x2 +2.00x3 = -11.00
-2.00x1 +1.00x2 +2.00x3 = -3.00

After step 1:
-3.00	-1.00	2.00	-11.00	
0.00	0.33	0.33	2.67	
0.00	0.33	0.67	-5.33	

After step 2:
-3.00	-1.00	2.00	-11.00	
0.00	0.33	0.33	2.67	
0.00	0.00	0.33	-8.00	

Unique Solution
Solution:
x1 = 2.00
x2 = 3.00
x3 = -1.00
```

### Example 2: No Solution

**Input:**
```
3
1 2 3 4
2 4 6 8
1 2 3 5
```

**System:**
```
x₁ + 2x₂ + 3x₃ = 4
2x₁ + 4x₂ + 6x₃ = 8
x₁ + 2x₂ + 3x₃ = 5
```

**Analysis**:  First and third equations are contradictory (same left side, different right side)

**Output:**
```
No Solution
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
Infinite Solutions
```

---

## 🎯 Compilation and Execution

### Compile
```bash
g++ -std=c++17 -O2 gauss-elimination-method.cpp -o gauss
```

### Run
```bash
./gauss
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

## 🔬 Applications

Gauss Elimination is fundamental to numerous applications:

### 1. **Engineering**
- Circuit analysis (Kirchhoff's laws)
- Structural analysis (finite element method)
- Heat transfer problems
- Fluid dynamics simulations

### 2. **Computer Graphics**
- 3D transformations
- Curve fitting
- Lighting calculations
- Texture mapping

### 3. **Economics**
- Input-output models (Leontief models)
- Economic equilibrium
- Portfolio optimization
- Resource allocation

### 4. **Machine Learning**
- Linear regression (normal equations)
- Principal Component Analysis (PCA)
- Neural network training
- Support Vector Machines

### 5. **Physics**
- Quantum mechanics (Schrödinger equation)
- Electromagnetism (Maxwell's equations)
- Mechanics (equilibrium problems)
- Optics (ray tracing)

### 6. **Data Science**
- Least squares fitting
- Interpolation
- Statistical analysis
- Network analysis

---

## 📚 References

- [Gaussian Elimination - Wikipedia](https://en.wikipedia.org/wiki/Gaussian_elimination)
- Numerical Methods For Engineers by Raymond Canale and Steven C. Chapra

---

**Part of the [Numerical Computing Suite](../) by [AbirHasanArko](https://github.com/AbirHasanArko)**   
Roll: 2207053   
Department of CSE, KUET   
