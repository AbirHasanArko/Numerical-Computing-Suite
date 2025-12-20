# Numerical Computing Suite

<p align="left">
  <img src="https://img.shields.io/badge/C++-Numerical%20Methods-blue?style=for-the-badge&logo=cplusplus" alt="C++ Numerical Methods"/>
  <a href="#authors">
    <img src="https://img.shields.io/badge/Authors-brown?style=for-the-badge&logo=github" alt="Authors"/>
  </a>
</p>


An extensive suite of C++ programs for solving a wide range of numerical computing problems, including linear and non-linear equations, interpolation, curve fitting, numerical differentiation, and integration. Each method is accompanied by theory, code, input/output samples, and usage notes.

---

## Table of Contents
  
- [Solution of Linear Equations](#solution-of-linear-equations)
    - [Gauss Elimination Method](#gauss-elimination-method)
        - [Theory](#gauss-elimination-theory)
        - [Code](#gauss-elimination-code)
        - [Input](#gauss-elimination-input)
        - [Output](#gauss-elimination-output)
    - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
        - [Theory](#gauss-jordan-theory)
        - [Code](#gauss-jordan-code)
        - [Input](#gauss-jordan-input)
        - [Output](#gauss-jordan-output)
    - [LU Decomposition](#lu-decomposition-method)
        - [Theory](#lu-decomposition-theory)
        - [Code](#lu-decomposition-code)
        - [Input](#lu-decomposition-input)
        - [Output](#lu-decomposition-output)
    - [Matrix Inversion](#matrix-inversion)
        - [Theory](#matrix-inversion-theory)
        - [Code](#matrix-inversion-code)
        - [Input](#matrix-inversion-input)
        - [Output](#matrix-inversion-output)

- [Solution of Non-linear Equations](#solution-of-non-linear-equations)
    - [Bisection Method](#bisection-method)
        - [Theory](#bisection-theory)
        - [Code](#bisection-code)
        - [Input](#bisection-input)
        - [Output](#bisection-output)
    - [False Position Method](#false-position-method)
        - [Theory](#false-position-theory)
        - [Code](#false-position-code)
        - [Input](#false-position-input)
        - [Output](#false-position-output)
    - [Newton-Raphson Method](#newton-raphson-method)
        - [Theory](#newton-raphson-theory)
        - [Code](#newton-raphson-code)
        - [Input](#newton-raphson-input)
        - [Output](#newton-raphson-output)
    - [Secant Method](#secant-method)
        - [Theory](#secant-theory)
        - [Code](#secant-code)
        - [Input](#secant-input)
        - [Output](#secant-output)

- [Interpolation and Approximation](#interpolation-and-approximation)
    - [Newton's Forward Interpolation](#newtons-forward-interpolation)
        - [Theory](#newtons-forward-theory)
        - [Code](#newtons-forward-code)
        - [Input](#newtons-forward-input)
        - [Output](#newtons-forward-output)
    - [Newton's Backward Interpolation](#newtons-backward-interpolation)
        - [Theory](#newtons-backward-theory)
        - [Code](#newtons-backward-code)
        - [Input](#newtons-backward-input)
        - [Output](#newtons-backward-output)
    - [Newton's Divided Difference Interpolation](#newtons-divided-difference-interpolation)
        - [Theory](#newtons-divided-difference-theory)
        - [Code](#newtons-divided-difference-code)
        - [Input](#newtons-divided-difference-input)
        - [Output](#newtons-divided-difference-output)

- [Numerical Integration](#numerical-integration)
    - [Simpson's One-third Rule](#simpsons-one-third-rule)
        - [Theory](#simpsons-one-third-theory)
        - [Code](#simpsons-one-third-code)
        - [Input](#simpsons-one-third-input)
        - [Output](#simpsons-one-third-output)
    - [Simpson's Three-eighths Rule](#simpsons-three-eighths-rule)
        - [Theory](#simpsons-three-eighths-theory)
        - [Code](#simpsons-three-eighths-code)
        - [Input](#simpsons-three-eighths-input)
        - [Output](#simpsons-three-eighths-output)

- [Numerical Differentiation](#numerical-differentiation)
    - [First and Second Order Derivative based on Forward Interpolation](#first-and-second-order-derivative-based-on-forward-interpolation)
        - [Theory](#forward-interpolation-derivative-theory)
        - [Code](#forward-interpolation-derivative-code)
        - [Input](#forward-interpolation-derivative-input)
        - [Output](#forward-interpolation-derivative-output)
    - [First and Second Order Derivative based on Backward Interpolation](#first-and-second-order-derivative-based-on-backward-interpolation)
        - [Theory](#backward-interpolation-derivative-theory)
        - [Code](#backward-interpolation-derivative-code)
        - [Input](#backward-interpolation-derivative-input)
        - [Output](#backward-interpolation-derivative-output)

- [Solution of Differential Equations](#solution-of-differential-equations)
    - [Runge Kutta](#runge-kutta)
        - [Theory](#runge-kutta-theory)
        - [Code](#runge-kutta-code)
        - [Input](#runge-kutta-input)
        - [Output](#runge-kutta-output)

- [Curve Fitting (Regression)](#curve-fitting-regression)
    - [Least Squares Lines](#least-squares-lines)
        - [Theory](#least-squares-lines-theory)
        - [Code](#least-squares-lines-code)
        - [Input](#least-squares-lines-input)
        - [Output](#least-squares-lines-output)
    - [Least Square Polynomials](#least-square-polynomials)
        - [Theory](#least-square-polynomials-theory)
        - [Code](#least-square-polynomials-code)
        - [Input](#least-square-polynomials-input)
        - [Output](#least-square-polynomials-output)
    - [Non-linear Curve Fitting](#non-linear-curve-fitting)
        - [Theory](#non-linear-curve-fitting-theory)
        - [Code](#non-linear-curve-fitting-code)
        - [Input](#non-linear-curve-fitting-input)
        - [Output](#non-linear-curve-fitting-output)
---
---

# Solution of Linear Equations
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-red?style=for-the-badge)](./Solution%20of%20Linear%20Equations/)
### 📐 Mathematical Foundation

#### What are Linear Equations?

A **linear equation** is an algebraic equation where each term is either a constant or the product of a constant and a single variable. The general form of a linear equation in *n* variables is:

```
a₁x₁ + a₂x₂ + a₃x₃ + ...  + aₙxₙ = b
```

Where:
- `a₁, a₂, ..., aₙ` are coefficients (constants)
- `x₁, x₂, ..., xₙ` are variables (unknowns)
- `b` is a constant term

**Properties of Linear Equations:**
- Each variable appears only to the first power
- Variables are not multiplied together (no x₁x₂ terms)
- No variables appear in denominators
- No transcendental functions (sin, cos, exp, log, etc.)

### Matrix Representation

A system of *n* linear equations with *n* unknowns can be represented in matrix form:

```
Ax = b
```

Where:
- **A** is an *n × n* coefficient matrix
- **x** is an *n × 1* solution vector (unknowns)
- **b** is an *n × 1* constant vector (right-hand side)

**Example**:  The system
```
2x + y - z = 8
-3x - y + 2z = -11
-2x + y + 2z = -3
```

Can be written as:
```
┌          ┐   ┌   ┐   ┌    ┐
│  2  1 -1 │   │ x │   │  8 │
│ -3 -1  2 │ × │ y │ = │-11 │
│ -2  1  2 │   │ z │   │ -3 │
└          ┘   └   ┘   └    ┘
```

**Augmented Matrix Representation**:  `[A|b]`
```
┌              ┐
│  2  1 -1 | 8 │
│ -3 -1  2 |-11│
│ -2  1  2 |-3 │
└              ┘
```

### Types of Solutions

A system of linear equations can have:  

#### 1. **Unique Solution** 🎯

**Definition**:  Exactly one set of values satisfies all equations simultaneously.

**Conditions**:
- The coefficient matrix **A** is non-singular (det(A) ≠ 0)
- rank(A) = rank([A|b]) = n
- All rows are linearly independent

**Geometric Interpretation**:
- In 2D: Two non-parallel lines intersect at one point
- In 3D:  Three non-parallel planes intersect at one point
- In nD:  Hyperplanes intersect at exactly one point

**Example**:
```
x + y = 3
x - y = 1
```
Solution: x = 2, y = 1 (lines intersect at (2, 1))

#### 2. **No Solution** ❌

**Definition**: No set of values satisfies all equations simultaneously. 

**Conditions**:
- The system is **inconsistent**
- rank(A) < rank([A|b])
- At least one equation contradicts the others

**Geometric Interpretation**:
- In 2D:  Parallel lines that never meet
- In 3D:  Planes that don't have a common intersection point
- System has contradictory constraints

**Example**:
```
x + y = 2
x + y = 5
```
No solution:  parallel lines (same slope, different intercepts)

#### 3. **Infinite Solutions** ∞

**Definition**:  Infinitely many sets of values satisfy the equations.

**Conditions**: 
- The system is **consistent but underdetermined**
- rank(A) = rank([A|b]) < n
- Some equations are linearly dependent (redundant)

**Geometric Interpretation**:
- In 2D: Coincident lines (same line)
- In 3D: Planes intersect along a line or coincide
- System has free variables (degrees of freedom)

**Example**:
```
x + y = 2
2x + 2y = 4
```
Infinite solutions: y = 2 - x (any point on the line)

### Fundamental Concepts

#### **Rank of a Matrix**

The **rank** is the maximum number of linearly independent rows (or columns) in a matrix. 

**Properties**:
- rank(A) ≤ min(m, n) for an m×n matrix
- rank(A) = n means A has full rank (non-singular for square matrices)
- rank(A) < n means A is singular (det(A) = 0)

**Calculation**:  After row reduction, count non-zero rows.

#### **Row Echelon Form (REF)**

A matrix is in row echelon form if:
1. All zero rows are at the bottom
2. The leading entry (pivot) of each non-zero row is to the right of the pivot in the row above
3. All entries below a pivot are zero

**Example**:
```
┌         ┐
│ 2  1 -1 │  ← pivot at position (1,1)
│ 0  3  2 │  ← pivot at position (2,2)
│ 0  0  5 │  ← pivot at position (3,3)
└         ┘
```

#### **Reduced Row Echelon Form (RREF)**

A matrix is in RREF if it's in REF and additionally:
1. All pivots are 1
2. Each pivot is the only non-zero entry in its column

**Example**:
```
┌         ┐
│ 1  0  0 │
│ 0  1  0 │
│ 0  0  1 │
└         ┘
```

#### **Determinant**

The determinant is a scalar value that encodes information about a square matrix. 

**Properties**:
- det(A) ≠ 0 ⟺ A is non-singular ⟺ unique solution exists
- det(A) = 0 ⟺ A is singular ⟺ infinite or no solutions
- det(AB) = det(A) × det(B)
- det(A⁻¹) = 1/det(A)   
   
---
   
# Gauss Elimination Method
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-green?style=for-the-badge)](./Solution%20of%20Linear%20Equations/Gauss-Jordan%20Elimination%20Method/)

## Gauss Elimination Theory
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
[  ⋮    ⋮     ⋱   ⋮   | ⋮  ]
[ aₙ₁  aₙ₂   ...  aₙₙ  | bₙ ]
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

## Gauss Elimination Code
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

---

## Gauss Elimination Input
### Input Format

```
n
a₁₁ a₁₂ ... a₁ₙ b₁
a₂₁ a₂₂ ...  a₂ₙ b₂
...
aₙ₁ aₙ₂ ... aₙₙ bₙ
```
**Input (input.txt):**   
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3

3
1 1 1 2
2 2 2 4
1 1 1 5

3
1 2 3 6
2 4 6 12
3 6 9 18
```

---

## Gauss Elimination Output
**Output (output.txt):** 
```
Input system:
2.00x1 +1.00x2 -1.00x3 = 8.00
-3.00x1 -1.00x2 +2.00x3 = -11.00
-2.00x1 +1.00x2 +2.00x3 = -3.00

After step 1:
-3.00	-1.00	2.00	-11.00	
0.00	0.33	0.33	0.67	
0.00	1.67	0.67	4.33	

After step 2:
-3.00	-1.00	2.00	-11.00	
0.00	1.67	0.67	4.33	
0.00	0.00	0.20	-0.20	

Unique Solution
Solution:
x1 = 2.00
x2 = 3.00
x3 = -1.00

-------------------------------------------


Input system:
1.00x1 +1.00x2 +1.00x3 = 2.00
2.00x1 +2.00x2 +2.00x3 = 4.00
1.00x1 +1.00x2 +1.00x3 = 5.00

After step 1:
2.00	2.00	2.00	4.00	
0.00	0.00	0.00	0.00	
0.00	0.00	0.00	3.00	

No Solution

-------------------------------------------


Input system:
1.00x1 +2.00x2 +3.00x3 = 6.00
2.00x1 +4.00x2 +6.00x3 = 12.00
3.00x1 +6.00x2 +9.00x3 = 18.00

After step 1:
3.00	6.00	9.00	18.00	
0.00	0.00	0.00	0.00	
0.00	0.00	0.00	0.00	

Infinite Solutions

-------------------------------------------
```

---

# Gauss-Jordan Elimination Method
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-blue?style=for-the-badge)](./Solution%20of%20Linear%20Equations/Gauss-Jordan%20Elimination%20Method/)

## Gauss Jordan Theory
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

## Gauss Jordan Code
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

## Gauss Jordan Input
### Input Format

```
n
a₁₁ a₁₂ ... a₁ₙ b₁
a₂₁ a₂₂ ...  a₂ₙ b₂
...
aₙ₁ aₙ₂ ... aₙₙ bₙ
```
**Input (input.txt):**   
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3

3
1 1 1 2
2 2 2 4
1 1 1 5

3
1 2 3 6
2 4 6 12
3 6 9 18

4
2 1 -1 1 3
1 3 2 -1 8
3 1 -3 2 5
1 2 1 -2 4
```

---

## Gauss Jordan Output
**Output (output.txt):** 
```

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

============================================================


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

============================================================


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

============================================================


Input system:
2.00x1 +1.00x2 -1.00x3 +1.00x4 = 3.00
1.00x1 +3.00x2 +2.00x3 -1.00x4 = 8.00
3.00x1 +1.00x2 -3.00x3 +2.00x4 = 5.00
1.00x1 +2.00x2 +1.00x3 -2.00x4 = 4.00

After swapping row 1 with row 3:
3.00	1.00	-3.00	2.00	5.00	
1.00	3.00	2.00	-1.00	8.00	
2.00	1.00	-1.00	1.00	3.00	
1.00	2.00	1.00	-2.00	4.00	

Step 1 - Making diagonal element a[1][1] = 1:
1.00	0.33	-1.00	0.67	1.67	
1.00	3.00	2.00	-1.00	8.00	
2.00	1.00	-1.00	1.00	3.00	
1.00	2.00	1.00	-2.00	4.00	

Step 2 - Eliminating column 1 in all other rows:
1.00	0.33	-1.00	0.67	1.67	
0.00	2.67	3.00	-1.67	6.33	
0.00	0.33	1.00	-0.33	-0.33	
0.00	1.67	2.00	-2.67	2.33	

Step 3 - Making diagonal element a[2][2] = 1:
1.00	0.33	-1.00	0.67	1.67	
0.00	1.00	1.12	-0.62	2.38	
0.00	0.33	1.00	-0.33	-0.33	
0.00	1.67	2.00	-2.67	2.33	

Step 4 - Eliminating column 2 in all other rows:
1.00	0.00	-1.38	0.88	0.88	
0.00	1.00	1.12	-0.62	2.38	
0.00	0.00	0.62	-0.12	-1.13	
0.00	0.00	0.12	-1.62	-1.63	

Step 5 - Making diagonal element a[3][3] = 1:
1.00	0.00	-1.38	0.88	0.88	
0.00	1.00	1.12	-0.62	2.38	
0.00	0.00	1.00	-0.20	-1.80	
0.00	0.00	0.12	-1.62	-1.63	

Step 6 - Eliminating column 3 in all other rows:
1.00	0.00	0.00	0.60	-1.60	
0.00	1.00	0.00	-0.40	4.40	
0.00	0.00	1.00	-0.20	-1.80	
0.00	0.00	0.00	-1.60	-1.40	

Step 7 - Making diagonal element a[4][4] = 1:
1.00	0.00	0.00	0.60	-1.60	
0.00	1.00	0.00	-0.40	4.40	
0.00	0.00	1.00	-0.20	-1.80	
0.00	0.00	0.00	1.00	0.88	

Step 8 - Eliminating column 4 in all other rows:
1.00	0.00	0.00	0.00	-2.13	
0.00	1.00	0.00	0.00	4.75	
0.00	0.00	1.00	0.00	-1.63	
0.00	0.00	0.00	1.00	0.88	

========================================
FINAL RESULT:
========================================

Unique Solution

Final Reduced Row Echelon Form (RREF):
1.00	0.00	0.00	0.00	-2.13	
0.00	1.00	0.00	0.00	4.75	
0.00	0.00	1.00	0.00	-1.63	
0.00	0.00	0.00	1.00	0.88	

Solution:
x1 = -2.13
x2 = 4.75
x3 = -1.63
x4 = 0.88

============================================================


```

---

# LU Decomposition Method
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Solution%20of%20Linear%20Equations/LU%20Decomposition/)

## LU Decomposition Theory
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

## LU Decomposition Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin)
    {
        cerr << "Error:  input.txt not found!\n";
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

        fout << "\nPerforming LU Decomposition...\n";

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
            // Continue decomposition as much as possible to check consistency
            
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
                
                if (!rowIsZero && fabs(U[i][i]) > 1e-12)
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

            // Verification: compute A*x and compare with b
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

## LU Decomposition Input
### Input Format

```
n
a₁₁ a₁₂ ... a₁ₙ b₁
a₂₁ a₂₂ ...  a₂ₙ b₂
...
aₙ₁ aₙ₂ ... aₙₙ bₙ
```
**Input (input.txt):**   
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3

3
1 1 1 2
2 2 2 4
1 1 1 5

3
1 2 3 6
2 4 6 12
3 6 9 18

4
2 1 -1 1 3
1 3 2 -1 8
3 1 -3 2 5
1 2 1 -2 4
```

---

## LU Decomposition Output
**Output (output.txt):** 
```

========================================
Input system:
2.0000x1 +1.0000x2 -1.0000x3 = 8.0000
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
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 
---------------------------------------------

After step 2:
L matrix:
    1.0000     0.0000     0.0000 
   -1.5000     1.0000     0.0000 
   -1.0000     4.0000     0.0000 
U matrix:
    2.0000     1.0000    -1.0000 
    0.0000     0.5000     0.5000 
    0.0000     0.0000     0.0000 
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


========================================
Input system:
2.0000x1 +1.0000x2 -1.0000x3 +1.0000x4 = 3.0000
1.0000x1 +3.0000x2 +2.0000x3 -1.0000x4 = 8.0000
3.0000x1 +1.0000x2 -3.0000x3 +2.0000x4 = 5.0000
1.0000x1 +2.0000x2 +1.0000x3 -2.0000x4 = 4.0000
========================================

Performing LU Decomposition...

After step 1:
L matrix:
    1.0000     0.0000     0.0000     0.0000 
    0.5000     0.0000     0.0000     0.0000 
    1.5000     0.0000     0.0000     0.0000 
    0.5000     0.0000     0.0000     0.0000 
U matrix:
    2.0000     1.0000    -1.0000     1.0000 
    0.0000     0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000     0.0000 
---------------------------------------------

After step 2:
L matrix:
    1.0000     0.0000     0.0000     0.0000 
    0.5000     1.0000     0.0000     0.0000 
    1.5000    -0.2000     0.0000     0.0000 
    0.5000     0.6000     0.0000     0.0000 
U matrix:
    2.0000     1.0000    -1.0000     1.0000 
    0.0000     2.5000     2.5000    -1.5000 
    0.0000     0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000     0.0000 
---------------------------------------------

After step 3:
L matrix:
    1.0000     0.0000     0.0000     0.0000 
    0.5000     1.0000     0.0000     0.0000 
    1.5000    -0.2000     1.0000     0.0000 
    0.5000     0.6000    -0.0000     0.0000 
U matrix:
    2.0000     1.0000    -1.0000     1.0000 
    0.0000     2.5000     2.5000    -1.5000 
    0.0000     0.0000    -1.0000     0.2000 
    0.0000     0.0000     0.0000     0.0000 
---------------------------------------------

After step 4:
L matrix:
    1.0000     0.0000     0.0000     0.0000 
    0.5000     1.0000     0.0000     0.0000 
    1.5000    -0.2000     1.0000     0.0000 
    0.5000     0.6000    -0.0000     1.0000 
U matrix:
    2.0000     1.0000    -1.0000     1.0000 
    0.0000     2.5000     2.5000    -1.5000 
    0.0000     0.0000    -1.0000     0.2000 
    0.0000     0.0000     0.0000    -1.6000 
---------------------------------------------

========================================
FINAL RESULT:
========================================

Unique Solution
Determinant of U = 8.0000

Final L matrix (Lower Triangular):
    1.0000     0.0000     0.0000     0.0000 
    0.5000     1.0000     0.0000     0.0000 
    1.5000    -0.2000     1.0000     0.0000 
    0.5000     0.6000    -0.0000     1.0000 

Final U matrix (Upper Triangular):
    2.0000     1.0000    -1.0000     1.0000 
    0.0000     2.5000     2.5000    -1.5000 
    0.0000     0.0000    -1.0000     0.2000 
    0.0000     0.0000     0.0000    -1.6000 

--- Forward Substitution (L*y = b) ---
y1 = 3.0000
y2 = 6.5000
y3 = 1.8000
y4 = -1.4000

--- Back Substitution (U*x = y) ---

Solution Vector (x):
x1 = -2.1250
x2 = 4.7500
x3 = -1.6250
x4 = 0.8750

--- Verification (A*x = b) ---
Row 1: 3.0000 (expected: 3.0000)
Row 2: 8.0000 (expected: 8.0000)
Row 3: 5.0000 (expected: 5.0000)
Row 4: 4.0000 (expected: 4.0000)

============================================================


```

---
# Matrix Inversion
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Solution%20of%20Linear%20Equations/Matrix%20Inversion/)

## Matrix Inversion Theory
### Mathematical Foundation

Given a square matrix $A$ of order $n$:

- The **adjugate** (adjoint) of $A$, denoted $\operatorname{adj}(A)$, is the transpose of its cofactor matrix.
- The **determinant** $\det(A)$ is a scalar value that determines if $A$ is invertible.
- The **inverse** of $A$ (if it exists) is:

$$
A^{-1} = \frac{1}{\det(A)} \operatorname{adj}(A)
$$

### Algorithm Steps

1. **Input**: Read $n$ and the $n \times n$ matrix $A$ from file.
2. **Determinant**: Compute $\det(A)$ using cofactor expansion.
3. **Check Singularity**: If $|\det(A)| < \varepsilon$ (very small), matrix is singular.
4. **Adjugate**: Compute $\operatorname{adj}(A)$ by calculating cofactors and transposing.
5. **Inverse**: Compute $A^{-1}$ as $\operatorname{adj}(A)/\det(A)$.
6. **Output**: Write the inverse matrix to file, or report if singular.

### Singularity Detection

A matrix is **singular** if $\det(A) = 0$. In practice, due to floating-point errors, a small threshold $\varepsilon$ (e.g., $10^{-12}$) is used:

- If $|\det(A)| < \varepsilon$, the matrix is considered singular and not invertible.

### Complexity Analysis

- **Determinant (cofactor expansion)**: $O(n!)$ (not efficient for large $n$)
- **Adjugate calculation**: $O(n^4)$ (each cofactor is a determinant)
- **Overall**: Suitable for small matrices (e.g., $n \leq 6$)
---

## Matrix Inversion Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double determinant(const vector<vector<double>>& A);
vector<vector<double>> adjoint(const vector<vector<double>>& A);

/* ---------------------------
   Compute inverse by 1/det(A) * adj(A)
----------------------------*/
bool inverseByAdjoint(const vector<vector<double>>& A,
                      vector<vector<double>>& inv) {
    int n = A.size();
    for (auto &r : A)
        if (r.size() != n) return false; // must be square

    double detA = determinant(A);
    if (fabs(detA) < 1e-12) return false; // non-invertible

    vector<vector<double>> adjA = adjoint(A);

    inv.assign(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adjA[i][j] / detA;

    return true;
}

/* ---------------------------
   Determinant (cofactor expansion)
----------------------------*/
double determinant(const vector<vector<double>>& A) {
    int n = A.size();
    if (n == 1) return A[0][0];
    if (n == 2)
        return A[0][0]*A[1][1] - A[0][1]*A[1][0];

    double det = 0;
    for (int col = 0; col < n; col++) {
        vector<vector<double>> sub(n-1, vector<double>(n-1));
        for (int i = 1; i < n; i++) {
            int c2 = 0;
            for (int j = 0; j < n; j++) {
                if (j == col) continue;
                sub[i-1][c2++] = A[i][j];
            }
        }
        det += ((col % 2 == 0) ? 1 : -1) * A[0][col] * determinant(sub);
    }
    return det;
}

/* ---------------------------
   Adjoint matrix
----------------------------*/
vector<vector<double>> adjoint(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> adj(n, vector<double>(n));

    if (n == 1) {
        adj[0][0] = 1;
        return adj;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vector<vector<double>> sub(n-1, vector<double>(n-1));
            int r2 = 0;
            for (int r = 0; r < n; r++) {
                if (r == i) continue;
                int c2 = 0;
                for (int c = 0; c < n; c++) {
                    if (c == j) continue;
                    sub[r2][c2++] = A[r][c];
                }
                r2++;
            }

            double cofactor = ((i + j) % 2 == 0 ? 1 : -1)
                              * determinant(sub);
            adj[j][i] = cofactor; // transpose
        }
    }
    return adj;
}

/* ---------------------------
   Main: File I/O
----------------------------*/
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin.is_open()) {
        cerr << "Failed to open input.txt\n";
        return 1;
    }

    int n;
    fin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fin >> A[i][j];

    vector<vector<double>> inv;
    fout << fixed << setprecision(6);

    if (!inverseByAdjoint(A, inv)) {
        fout << "Matrix is singular.\n";
        return 0;
    }

    fout << "Inverse using adj(A)/det(A):\n";
    for (auto &row : inv) {
        for (double x : row)
            fout << setw(12) << x;
        fout << "\n";
    }

    fin.close();
    fout.close();
    return 0;
}
```

---

## Matrix Inversion Input
**Input (input.txt):**   
```
3
4 7 2
3 6 1
2 5 3
```
---

## Matrix Inversion Output
**Output (output.txt):** 
```
Inverse using adj(A)/det(A):
    1.444444   -1.222222   -0.555556
   -0.777778    0.888889    0.222222
    0.333333   -0.666667    0.333333

```
---
---

# Solution of Non-Linear Equations
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-red?style=for-the-badge)](./Solution%20of%20Linear%20Equations/)
## 📐 Mathematical Foundation

### What is a Non-linear Equation?

A **non-linear equation** is any equation where the variable $x$ appears with exponents other than one, inside transcendental functions (sin, exp, log, etc.), or in products with itself. Examples:
- $x^3 - x - 2 = 0$
- $e^x + x = 5$
- $\sin(x) = x/2$

### Root-finding Problem

The goal is to find $x^*$ such that $f(x^*) = 0$. Analytical solutions are rare, so we use iterative numerical methods:
- **Bracketing methods**: Start with an interval $[a, b]$ where $f(a)$ and $f(b)$ have opposite signs (guaranteed root by Intermediate Value Theorem).
- **Open methods**: Use one or two initial guesses and iterate using function values (and possibly derivatives).

### Convergence and Error

- **Convergence**: How quickly the method approaches the root.
- **Order of convergence**: Linear (slow), superlinear, quadratic (fast).
- **Stopping criteria**: $|f(x_n)| < \epsilon$ or $|x_{n+1} - x_n| < \epsilon$ for some small $\epsilon$.
- **Error estimation**: Each method provides a way to estimate or bound the error at each step.
---
   
# Bisection Method
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-cyan?style=for-the-badge)](./Solution%20of%20Non-linear%20Equations/Bisection%20Method/)

## Bisection Theory
### Mathematical Foundation

The Bisection Method is based on the **Intermediate Value Theorem** (IVT), which states:

> If a continuous function $f(x)$ satisfies $f(a) \cdot f(b) < 0$ on an interval $[a, b]$, then there exists at least one root $r \in (a, b)$ such that $f(r) = 0$.

**Core Principle:**

Given an interval $[x_{\text{low}}, x_{\text{high}}]$ where $f(x_{\text{low}}) \cdot f(x_{\text{high}}) < 0$:

1. Calculate the midpoint: $x_{\text{mid}} = \frac{x_{\text{low}} + x_{\text{high}}}{2}$

2. Evaluate $f(x_{\text{mid}})$:
   - If $f(x_{\text{mid}}) \approx 0$: Root found at $x_{\text{mid}}$
   - If $f(x_{\text{low}}) \cdot f(x_{\text{mid}}) < 0$: Root is in $[x_{\text{low}}, x_{\text{mid}}]$, set $x_{\text{high}} = x_{\text{mid}}$
   - If $f(x_{\text{mid}}) \cdot f(x_{\text{high}}) < 0$: Root is in $[x_{\text{mid}}, x_{\text{high}}]$, set $x_{\text{low}} = x_{\text{mid}}$

3. Repeat until convergence criteria are met

**Convergence Criteria:**

The algorithm stops when either condition is satisfied:

$$|f(x_{\text{mid}})| < \epsilon_f \quad \text{or} \quad |x_{\text{high}} - x_{\text{low}}| < \epsilon_x$$

where:
- $\epsilon_f = 10^{-6}$ (function value tolerance)
- $\epsilon_x = 10^{-6}$ (interval width tolerance)

**Error Bound:**

After $n$ iterations, the error is bounded by:

$$|r - x_{\text{mid}}| \leq \frac{b - a}{2^{n+1}}$$

where $[a, b]$ is the initial interval containing the root.

### Algorithm Steps

1. **Input Processing**:
   - Read polynomial degree $n$
   - Read $n+1$ coefficients: $a_0, a_1, ..., a_n$ for $f(x) = a_0 x^n + a_1 x^{n-1} + ... + a_{n-1}x + a_n$

2. **Interval Scanning**:
   - Define search range: $[-R, R]$ where $R = 5000$ (default)
   - Use step size: $h = 0.5$
   - For each point $x_i = -R + i \cdot h$:
     - Check if $|f(x_i)| < \epsilon$ (direct root detection)
     - Check if $f(x_i) \cdot f(x_{i+1}) < 0$ (sign change detection)
     - Store intervals with sign changes

3. **Bisection Iteration** (for each interval):
   - Initialize: $x_{\text{low}} = x_i$, $x_{\text{high}} = x_{i+1}$
   - While not converged:
     1. $x_{\text{mid}} = \frac{x_{\text{low}} + x_{\text{high}}}{2}$
     2. $f_{\text{mid}} = f(x_{\text{mid}})$
     3. Record iteration: $(iteration, x_{\text{low}}, x_{\text{high}}, x_{\text{mid}}, f_{\text{mid}})$
     4. Check convergence: $|f_{\text{mid}}| < \epsilon$ or $|x_{\text{high}} - x_{\text{low}}| < \epsilon$
     5. Update interval:
        - If $f(x_{\text{low}}) \cdot f_{\text{mid}} < 0$: $x_{\text{high}} = x_{\text{mid}}$
        - Else: $x_{\text{low}} = x_{\text{mid}}$

4. **Root Collection**:
   - Add converged $x_{\text{mid}}$ to roots list
   - Check for duplicates (within tolerance $\epsilon$)
   - Sort roots in ascending order

5. **Output Generation**:
   - Display polynomial equation
   - Print final roots table
   - Print iteration history for each root

### Complexity Analysis

- **Time Complexity**:
  - Interval scanning: $O(R/h)$ where $R$ is search range, $h$ is step size
  - Single root finding: $O(\log_2(\frac{b-a}{\epsilon}))$ where $[a,b]$ is initial interval
  - For $m$ roots: $O(R/h + m \log_2(\frac{b-a}{\epsilon}))$
  - Typical: $O(R/h)$ dominates for large search ranges

- **Space Complexity**: $O(m \cdot k)$ where $m$ is number of roots, $k$ is average iterations per root

- **Convergence Rate**: Linear convergence - error reduces by half each iteration

- **Accuracy**: 
  - Guaranteed to converge to a root within tolerance $\epsilon$
  - Final error: $|r_{\text{true}} - r_{\text{computed}}| \lesssim \epsilon$

---

## Bisection Code
```cpp
#include <bits/stdc++.h>
using namespace std;

#define f double

// Polynomial function f(x)
f fun(const vector<f> &coef, f x) {
    f result = 0.0;
    int n = coef.size() - 1;
    for (int i = 0; i <= n; i++)
        result += coef[i] * pow(x, n - i);
    return result;
}

// Print polynomial nicely
void printPolynomial(ofstream &out, const vector<f> &coef) {
    int n = coef.size() - 1;
    out << "f(x) = ";
    bool first = true;
    for (int i = 0; i <= n; i++) {
        f c = coef[i]; int power = n - i;
        if (c == 0) continue;
        if (!first) {
            if (c > 0) out << " + "; else out << " - ";
        } else { if (c < 0) out << "-"; }
        if (abs(c) != 1 || power == 0) out << abs(c);
        if (power > 0) { out << "x"; if (power > 1) out << "^" << power; }
        first = false;
    }
    out << "\n\n";
}

// Print header
void printHeader(ofstream &out) {
    out << "============================================\n";
    out << "              Bisection Method\n";
    out << "============================================\n\n";
}
void printHeader() {
    cout << "============================================\n";
    cout << "              Bisection Method\n";
    cout << "============================================\n\n";
}

// Print final roots table with trivial root highlight
void printRootsTable(ofstream &out, const vector<f> &roots, const vector<f> &coef, f tolerance = 1e-6) {
    out << setw(6) << "Index" 
        << setw(20) << "Root Value" << "\n";
    out << string(26, '-') << "\n";

    out << fixed << setprecision(6);
    for (size_t i = 0; i < roots.size(); i++) {
        f fx = fun(coef, roots[i]);
        out << setw(6) << (i + 1)
            << setw(20) << roots[i]
            << "\n";
    }
    out << "\n";
}

// Print iteration table per root
void printBisectionIteration(ofstream &out, const vector<tuple<int,f,f,f,f>> &iterData) {
    out << setw(10) << "Iteration"
        << setw(15) << "x_low"
        << setw(15) << "x_high"
        << setw(15) << "x_mid"
        << setw(20) << "f(x_mid)" << "\n";
    out << string(75, '-') << "\n";
    out << fixed << setprecision(6);
    for (auto &t : iterData) {
        int iter; f xL, xH, xM, fxM;
        tie(iter, xL, xH, xM, fxM) = t;
        out << setw(10) << iter
            << setw(15) << xL
            << setw(15) << xH
            << setw(15) << xM
            << setw(20) << fxM << "\n";
    }
    out << "\n";
}

int main() {
    string inputFile, outputFile;

    // File names
    printHeader();
    cout << "Enter input file name: ";
    cin >> inputFile;
    cout << "Enter output file name: ";
    cin >> outputFile;

    ifstream in(inputFile);
    ofstream out(outputFile);
    if (!in.is_open() || !out.is_open()) {
        cout << "Error: Cannot open input/output file!\n";
        return 1;
    }

    int degree;
    in >> degree;
    vector<f> coef(degree + 1);
    for (int i = 0; i <= degree; i++) in >> coef[i];

    printHeader(out);
    out << "Polynomial Degree : " << degree << "\n";
    out << "Coefficients      : ";
    for (auto c : coef) out << c << " ";
    out << "\n\n";
    printPolynomial(out, coef);

    // Bisection parameters
    f searchRange = 5000.0, step = 0.5, tolerance = 1e-6;
    vector<f> roots, intervals;

    // Scan for roots / intervals
    for (f i = -searchRange; i <= searchRange; i += step) {
        f fx1 = fun(coef, i);
        f fx2 = fun(coef, i + step);
        if (abs(fx1) < tolerance) {
            bool duplicate = false;
            for (auto r : roots) if (abs(r - i) < tolerance) { duplicate = true; break; }
            if (!duplicate) roots.push_back(i);
        } else if (fx1 * fx2 < 0.0) intervals.push_back(i);
    }

    vector<vector<tuple<int,f,f,f,f>>> allIterations; // iterations per root

    for (f start : intervals) {
        f xL = start, xH = start + step, xM;
        vector<tuple<int,f,f,f,f>> iterData;
        int iter = 1;
        while (true) {
            xM = (xL + xH) / 2.0;
            f fxM = fun(coef, xM);
            iterData.push_back({iter, xL, xH, xM, fxM});

            if (abs(fxM) < tolerance || abs(xH - xL) < tolerance) {
                bool duplicate = false;
                for (auto r : roots) if (abs(r - xM) < tolerance) { duplicate = true; break; }
                if (!duplicate) roots.push_back(xM);
                break;
            }

            f fxL = fun(coef, xL);
            if (fxL * fxM < 0) xH = xM;
            else xL = xM;

            iter++;
        }
        allIterations.push_back(iterData);
    }

    sort(roots.begin(), roots.end());

    // Print final roots
    out << "Roots\n";
    out << "============================================\n";
    if (roots.empty()) {
        out << "No real roots found in the given range.\n\n";
    } else {
        printRootsTable(out, roots, coef, tolerance);
    }

    // Print iteration tables
    for (size_t i = 0; i < allIterations.size(); i++) {
        out << "Iteration Table for Root " << (i + 1) << "\n";
        out << "----------------------------------------\n";
        printBisectionIteration(out, allIterations[i]);
    }

    out << "============================================\n";
    out << "Computation Completed Successfully.\n";

    in.close();
    out.close();
    cout << "Computation completed. Results written to '" << outputFile << "'\n";

    return 0;
}
```

---

## Bisection Input
**Input1 (input1.txt):**   
```
3
1 -6 11 -6
```
**Input2 (input2.txt):**   
```
2
1 0 1
```
---

## Bisection Output
**Output1 (output1.txt):** 
```
============================================
              Bisection Method
============================================

Polynomial Degree : 3
Coefficients      : 1 -6 11 -6 

f(x) = x^3 - 6x^2 + 11x - 6

Roots
============================================
 Index          Root Value
--------------------------
     1            1.000000
     2            2.000000
     3            3.000000

============================================
Computation Completed Successfully.
```

**Output2 (output2.txt):** 
```
============================================
              Bisection Method
============================================

Polynomial Degree : 2
Coefficients      : 1 0 1 

f(x) = x^2 + 1

Roots
============================================
No real roots found in the given range.

============================================
Computation Completed Successfully.
```

---

# False Position Method
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-blue?style=for-the-badge)](./Solution%20of%20Non-linear%20Equations/False-Position%20Method/)

## False Position Theory
### Mathematical Foundation

The False Position Method is based on the **Intermediate Value Theorem** and uses **linear interpolation** to find the root. Unlike the Bisection Method which simply takes the midpoint, the False Position Method uses a weighted average that considers the function values.

**Core Principle:**

Given an interval $[x_L, x_R]$ where $f(x_L) \cdot f(x_R) < 0$, the method approximates the root by drawing a straight line (secant) connecting the points $(x_L, f(x_L))$ and $(x_R, f(x_R))$, and finding where this line crosses the x-axis.

**False Position Formula:**

$$x_0 = \frac{x_L \cdot f(x_R) - x_R \cdot f(x_L)}{f(x_R) - f(x_L)}$$

Alternatively, this can be written as:

$$x_0 = x_R - f(x_R) \cdot \frac{x_R - x_L}{f(x_R) - f(x_L)}$$

or as a weighted average:

$$x_0 = \frac{x_L \cdot |f(x_R)| + x_R \cdot |f(x_L)|}{|f(x_R)| + |f(x_L)|}$$

**Geometric Interpretation:**

The formula finds the x-intercept of the line passing through $(x_L, f(x_L))$ and $(x_R, f(x_R))$:

- If $|f(x_L)| > |f(x_R)|$: The root is closer to $x_R$, so $x_0$ is weighted toward $x_R$
- If $|f(x_L)| < |f(x_R)|$: The root is closer to $x_L$, so $x_0$ is weighted toward $x_L$
- If $|f(x_L)| = |f(x_R)|$: Reduces to midpoint (like bisection)

**Iteration Process:**

1. Calculate $x_0$ using the False Position formula
2. Evaluate $f(x_0)$
3. Update the interval:
   - If $f(x_L) \cdot f(x_0) < 0$: Root is in $[x_L, x_0]$, set $x_R = x_0$, $f(x_R) = f(x_0)$
   - If $f(x_0) \cdot f(x_R) < 0$: Root is in $[x_0, x_R]$, set $x_L = x_0$, $f(x_L) = f(x_0)$
4. Repeat until convergence

**Convergence Criteria:**

The algorithm stops when either condition is satisfied:

$$|f(x_0)| < \epsilon_f \quad \text{or} \quad |x_R - x_L| < \epsilon_x$$

where:
- $\epsilon_f = 10^{-6}$ (function value tolerance)
- $\epsilon_x = 10^{-6}$ (interval width tolerance)

**Comparison with Bisection:**

| Aspect | Bisection Method | False Position Method |
|--------|------------------|----------------------|
| Point Selection | $x_{\text{mid}} = \frac{x_L + x_R}{2}$ | $x_0 = \frac{x_L f(x_R) - x_R f(x_L)}{f(x_R) - f(x_L)}$ |
| Uses Function Values | No (only signs) | Yes (weighted by values) |
| Convergence Rate | Linear (guaranteed halving) | Superlinear (typically faster) |
| Worst Case | Predictable | Can be slow for certain functions |
| Best For | Any continuous function | Nearly linear functions |

### Algorithm Steps

1. **Input Processing**:
   - Read polynomial degree $n$
   - Read $n+1$ coefficients: $a_0, a_1, ..., a_n$ for $f(x) = a_0 x^n + a_1 x^{n-1} + ... + a_{n-1}x + a_n$

2. **Interval Scanning**:
   - Define search range: $[-R, R]$ where $R = 5000$ (default)
   - Use step size: $h = 0.5$
   - For each point $x_i = -R + i \cdot h$:
     - Check if $|f(x_i)| < \epsilon$ (direct root detection)
     - Check if $f(x_i) \cdot f(x_{i+1}) < 0$ (sign change detection)
     - Store intervals with sign changes

3. **False Position Iteration** (for each interval):
   - Initialize: $x_L = x_i$, $x_R = x_{i+1}$
   - Calculate: $f_L = f(x_L)$, $f_R = f(x_R)$
   - Verify: $f_L \cdot f_R < 0$ (must have opposite signs)
   - While not converged:
     1. Calculate new estimate: $x_0 = \frac{x_L \cdot f_R - x_R \cdot f_L}{f_R - f_L}$
     2. Evaluate: $f_0 = f(x_0)$
     3. Check convergence: $|f_0| < \epsilon$ or $|x_R - x_L| < \epsilon$
     4. If converged:
        - Check for duplicate in roots list
        - Add $x_0$ to roots if unique
        - Break loop
     5. Update interval:
        - If $f_L \cdot f_0 < 0$: Set $x_R = x_0$, $f_R = f_0$ (root is in left half)
        - Else: Set $x_L = x_0$, $f_L = f_0$ (root is in right half)

4. **Root Collection**:
   - Combine roots from direct detection and False Position iterations
   - Filter duplicates (within tolerance $\epsilon$)
   - Sort roots in ascending order

5. **Output Generation**:
   - Display polynomial equation
   - Print final roots table with indices
   - Report completion status

### Complexity Analysis

- **Time Complexity**:
  - Interval scanning: $O(R/h)$ where $R$ is search range, $h$ is step size
  - Single root finding: $O(k)$ where $k$ is number of iterations (typically less than bisection)
  - For $m$ roots: $O(R/h + m \cdot k)$
  - Average case: Faster than bisection due to better point selection

- **Space Complexity**: $O(m)$ where $m$ is number of roots (no iteration history stored)

- **Convergence Rate**: 
  - **Superlinear** for well-behaved functions
  - Can approach **quadratic** convergence when function is nearly linear near root
  - Typically requires **30-50% fewer iterations** than bisection

- **Accuracy**: 
  - Guaranteed to converge to a root within tolerance $\epsilon$
  - Final error: $|r_{\text{true}} - r_{\text{computed}}| \lesssim \epsilon$
  - Often achieves higher accuracy faster than bisection

---

## False Position Code
```cpp
#include <bits/stdc++.h>
using namespace std;

#define f double

// Polynomial function f(x)
f fun(const vector<f> &coef, f x) {
    f result = 0.0;
    int n = coef.size() - 1;
    for (int i = 0; i <= n; i++)
        result += coef[i] * pow(x, n - i);
    return result;
}

// Print polynomial in readable form
void printPolynomial(ofstream &out, const vector<f> &coef) {
    int n = coef.size() - 1;
    out << "f(x) = ";

    bool first = true;
    for (int i = 0; i <= n; i++) {
        f c = coef[i];
        int power = n - i;
        if (c == 0) continue;

        if (!first) {
            if (c > 0) out << " + ";
            else out << " - ";
        } else {
            if (c < 0) out << "-";
        }

        if (abs(c) != 1 || power == 0)
            out << abs(c);

        if (power > 0) {
            out << "x";
            if (power > 1)
                out << "^" << power;
        }

        first = false;
    }
    out << "\n\n";
}

// Print header
void printHeader(ofstream &out) {
    out << "============================================\n";
    out << "           False Position Method\n";
    out << "============================================\n\n";
}
void printHeader() {
    cout << "============================================\n";
    cout << "           False Position Method\n";
    cout << "============================================\n\n";
}

// Print roots table with trivial root highlight
void printRootsTable(ofstream &out, const vector<f> &roots, const vector<f> &coef, f tolerance = 1e-6) {
    out << setw(6) << "Index" 
        << setw(20) << "Root Value" << "\n";
    out << string(26, '-') << "\n";

    out << fixed << setprecision(6);
    for (size_t i = 0; i < roots.size(); i++) {
        f fx = fun(coef, roots[i]);
        out << setw(6) << (i + 1)
            << setw(20) << roots[i]
            << "\n";
    }
    out << "\n";
}

int main() {
    string inputFile, outputFile;

    printHeader();
    // File names
    cout << "Enter input file name: ";
    cin >> inputFile;
    cout << "Enter output file name: ";
    cin >> outputFile;

    ifstream in(inputFile);
    ofstream out(outputFile);

    if (!in.is_open() || !out.is_open()) {
        cout << "Error: Cannot open input/output file!\n";
        return 1;
    }

    int degree;
    in >> degree;
    vector<f> coef(degree + 1);
    for (int i = 0; i <= degree; i++)
        in >> coef[i];

    printHeader(out);

    out << "Polynomial Degree : " << degree << "\n";
    out << "Coefficients      : ";
    for (auto c : coef) out << c << " ";
    out << "\n\n";

    printPolynomial(out, coef);

    // Parameters
    f searchRange = 5000.0, step = 0.5, tolerance = 1e-6;
    vector<f> roots, intervals;

    // Find intervals with sign changes
    for (f i = -searchRange; i <= searchRange; i += step) {
        f fx1 = fun(coef, i);
        f fx2 = fun(coef, i + step);

        // Exact root at grid
        if (abs(fx1) < tolerance) {
            bool duplicate = false;
            for (auto r : roots) if (abs(r - i) < tolerance) { duplicate = true; break; }
            if (!duplicate) roots.push_back(i);
        }
        // Sign change
        else if (fx1 * fx2 < 0.0) {
            intervals.push_back(i);
        }
    }

    // Apply False Position Method
    for (f start : intervals) {
        f xL = start, xR = start + step, x0;
        f fL = fun(coef, xL), fR = fun(coef, xR);

        // Make sure fL*fR < 0
        if (fL * fR > 0) continue;

        while (true) {
            // False Position formula
            x0 = (xL*fR - xR*fL) / (fR - fL);
            f f0 = fun(coef, x0);

            if (abs(f0) < tolerance || abs(xR - xL) < tolerance) {
                bool duplicate = false;
                for (auto r : roots) if (abs(r - x0) < tolerance) { duplicate = true; break; }
                if (!duplicate) roots.push_back(x0);
                break;
            }

            if (fL * f0 < 0) {
                xR = x0;
                fR = f0;
            } else {
                xL = x0;
                fL = f0;
            }
        }
    }

    sort(roots.begin(), roots.end());

    // Output roots table
    out << "Roots\n";
    out << "============================================\n";
    if (roots.empty()) {
        out << "No real roots found in the given range.\n\n";
    } else {
        printRootsTable(out, roots, coef, tolerance);
    }

    out << "============================================\n";
    out << "Computation Completed Successfully.\n";

    in.close();
    out.close();

    cout << "Computation completed. Results written to '" << outputFile << "'\n";

    return 0;
}

```

---

## False Position Input
**Input1 (input1.txt):**   
```
3
1 -6 11 -6
```

**Input2 (input2.txt):**   
```
2
1 0 1
```

---

## False Position Output
**Output1 (output1.txt):** 
```
============================================
           False Position Method
============================================

Polynomial Degree : 3
Coefficients      : 1 -6 11 -6 

f(x) = x^3 - 6x^2 + 11x - 6

Roots
============================================
 Index          Root Value
--------------------------
     1            1.000000
     2            2.000000
     3            3.000000

============================================
Computation Completed Successfully.
```
**Output2 (output2.txt):** 
```
============================================
           False Position Method
============================================

Polynomial Degree : 2
Coefficients      : 1 0 1 

f(x) = x^2 + 1

Roots
============================================
No real roots found in the given range.

============================================
Computation Completed Successfully.
```
---

# Newton-Raphson Method
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Solution%20of%20Non-linear%20Equations/Newton-Raphson%20Method/)

## Newton-Raphson Theory
### Mathematical Foundation

The Newton-Raphson Method is based on **linear approximation** using the tangent line. It uses the first-order Taylor series expansion to approximate the root.

**Core Principle:**

Given a differentiable function $f(x)$ and an initial guess $x_0$, the method approximates the root by finding where the tangent line at $(x_0, f(x_0))$ crosses the x-axis.

**Newton-Raphson Formula:**

$$x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$$

where:
- $x_n$ is the current approximation
- $f(x_n)$ is the function value at $x_n$
- $f'(x_n)$ is the derivative (slope of tangent) at $x_n$
- $x_{n+1}$ is the next (improved) approximation

**Geometric Interpretation:**

At point $(x_n, f(x_n))$, the tangent line has equation:

$$y - f(x_n) = f'(x_n)(x - x_n)$$

Setting $y = 0$ to find the x-intercept:

$$0 - f(x_n) = f'(x_n)(x - x_n)$$

$$x = x_n - \frac{f(x_n)}{f'(x_n)}$$

This x-intercept becomes the next approximation $x_{n+1}$.

**Derivation from Taylor Series:**

The Taylor series expansion of $f(x)$ around $x_n$ is:

$$f(x) = f(x_n) + f'(x_n)(x - x_n) + \frac{f''(x_n)}{2}(x - x_n)^2 + ...$$

For the root $r$ where $f(r) = 0$, ignoring higher-order terms:

$$0 \approx f(x_n) + f'(x_n)(r - x_n)$$

$$r \approx x_n - \frac{f(x_n)}{f'(x_n)}$$

This gives the Newton-Raphson iteration formula.

**Iteration Process:**

1. Start with an initial guess $x_0$ (preferably near the root)
2. Compute function value: $f(x_n)$
3. Compute derivative: $f'(x_n)$
4. Update approximation: $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$
5. Check convergence: $|f(x_{n+1})| < \epsilon$
6. Repeat from step 2 with $x_n = x_{n+1}$ until convergence

**Convergence Criteria:**

The algorithm stops when:

$$|f(x_{n+1})| < \epsilon$$

where:
- $\epsilon = 10^{-6}$ (function value tolerance)

Additionally, the implementation includes **cycle detection** using a map to prevent infinite loops when the method oscillates.

**Convergence Rate - Quadratic:**

The Newton-Raphson Method has **quadratic convergence** when near the root, meaning:

$$|e_{n+1}| \approx C |e_n|^2$$

where $e_n = x_n - r$ is the error at iteration $n$. This means:
- **Error squares** each iteration
- **Number of correct digits doubles** each iteration
- **Extremely fast convergence** near the root

**Example of Quadratic Convergence:**
- Iteration 1: Error = $10^{-1}$ (1 digit)
- Iteration 2: Error = $10^{-2}$ (2 digits)
- Iteration 3: Error = $10^{-4}$ (4 digits)
- Iteration 4: Error = $10^{-8}$ (8 digits)
- Iteration 5: Error = $10^{-16}$ (machine precision)

**Comparison with Other Methods:**

| Aspect | Bisection | False Position | Secant | Newton-Raphson |
|--------|-----------|----------------|--------|----------------|
| Derivative Required | No | No | No | **Yes** |
| Convergence Order | 1.0 | ~1.5 | 1.618 | **2.0** |
| Iterations (typical) | ~20 | ~8-10 | ~4-6 | **~3-4** |
| Function Evals/Iter | 1 | 2 | 2 | 1 + derivative |
| Initial Points | 2 (bracket) | 2 (bracket) | 2 (any) | **1** |
| Bracketing Required | Yes | Yes | No | **No** |
| Speed | Slow | Medium | Fast | **Fastest** |
| Reliability | Guaranteed | High | Medium | **Medium** |
| Best For | Guaranteed | Linear funcs | No derivative | **Speed + derivative** |

**When Newton-Raphson Excels:**
- ✅ Derivative is easily computable (polynomials, exponentials, trig functions)
- ✅ High precision needed quickly
- ✅ Function is smooth and well-behaved near root
- ✅ Good initial guess available
- ✅ Computational resources allow derivative calculation

**When Newton-Raphson May Fail:**
- ❌ Derivative is zero or very small (horizontal tangent) → division by zero
- ❌ Poor initial guess far from root → may diverge or oscillate
- ❌ Multiple roots or inflection points nearby → may converge to wrong root
- ❌ Function has discontinuities or sharp corners → derivative undefined

### Algorithm Steps

1. **Input Processing**:
   - Read polynomial degree $n$
   - Read $n+1$ coefficients: $a_0, a_1, ..., a_n$ for $f(x) = a_0 x^n + a_1 x^{n-1} + ... + a_{n-1}x + a_n$

2. **Derivative Setup**:
   - Compute derivative coefficients using power rule
   - For polynomial $f(x) = \sum_{i=0}^{n} a_i x^{n-i}$
   - Derivative: $f'(x) = \sum_{i=0}^{n-1} a_i (n-i) x^{n-i-1}$

3. **Interval Scanning**:
   - Define search range: $[-R, R]$ where $R = 5000$ (default)
   - Use step size: $h = 0.5$
   - For each point $x_i = -R + i \cdot h$:
     - Check if $|f(x_i)| < \epsilon$ (direct root detection)
     - Check if $f(x_i) \cdot f(x_{i+1}) < 0$ (sign change detection)
     - Store intervals with sign changes as starting points

4. **Newton-Raphson Iteration** (for each starting point):
   - Initialize: $x = $ starting point
   - Create empty visited map for cycle detection
   - While not converged:
     1. Evaluate: $f(x)$ and $f'(x)$
     2. **Check for zero derivative:** If $|f'(x)| < 10^{-10}$, break (avoid division by zero)
     3. **Apply Newton-Raphson formula:** 
        $$x_{\text{new}} = x - \frac{f(x)}{f'(x)}$$
     4. Evaluate: $f(x_{\text{new}})$
     5. Round $x_{\text{new}}$ for map: $x_{\text{new,rounded}} = \text{round}(x_{\text{new}} \times 10^6) / 10^6$
     6. **Convergence/cycle check:**
        - If $|f(x_{\text{new}})| < \epsilon$ OR $x_{\text{new,rounded}}$ is in visited map:
          - Check for duplicate in roots list
          - Add $x_{\text{new}}$ to roots if unique
          - Break loop
     7. Mark $x_{\text{new,rounded}}$ as visited
     8. **Update for next iteration:** $x = x_{\text{new}}$

5. **Root Collection**:
   - Combine roots from direct detection and Newton-Raphson iterations
   - Filter duplicates (within tolerance $\epsilon$)
   - Sort roots in ascending order

6. **Output Generation**:
   - Display polynomial equation
   - Print final roots table with indices
   - Report completion status

### Complexity Analysis

- **Time Complexity**:
  - Interval scanning: $O(R/h)$ where $R$ is search range, $h$ is step size
  - Single root finding: $O(k)$ where $k$ is number of iterations (typically 3-5)
  - For $m$ roots: $O(R/h + m \cdot k)$
  - **Fastest convergence** among all methods: typically requires **60-80% fewer iterations** than Bisection

- **Space Complexity**: 
  - $O(m + k)$ where $m$ is number of roots
  - $k$ is average size of visited map per root (cycle detection overhead)
  - Negligible for most practical problems

- **Convergence Rate**: 
  - **Quadratic** with order 2.0 (when near root and $f''(r) \neq 0$)
  - Error reduces by factor of approximately $e^2$ each iteration
  - **Significantly faster** than all other methods:
    - ~6x faster than Bisection
    - ~3x faster than False Position
    - ~2x faster than Secant

- **Function Evaluations per Iteration**:
  - **1 function evaluation** $f(x)$
  - **1 derivative evaluation** $f'(x)$
  - If derivative is cheap (polynomials), very efficient
  - If derivative is expensive (complex functions), may be costlier than Secant

- **Accuracy**: 
  - Achieves machine precision ($\sim 10^{-16}$ for double) in just 5-6 iterations
  - Final error: $|r_{\text{true}} - r_{\text{computed}}| \approx 10^{-15}$ (limited by floating-point precision)
  - **Highest accuracy** in fewest iterations

---

## Newton-Raphson Code
```cpp
#include <bits/stdc++.h>
using namespace std;

#define f double

// Polynomial function f(x)
f fun(const vector<f> &coef, f x) {
    f result = 0.0;
    int n = coef.size() - 1;
    for (int i = 0; i <= n; i++)
        result += coef[i] * pow(x, n - i);
    return result;
}

// Derivative function f'(x)
f derivative(const vector<f> &coef, f x) {
    f result = 0.0;
    int n = coef.size() - 1;
    for (int i = 0; i < n; i++)
        result += coef[i] * (n - i) * pow(x, n - i - 1);
    return result;
}

// Print polynomial in readable form
void printPolynomial(ofstream &out, const vector<f> &coef) {
    int n = coef.size() - 1;
    out << "f(x) = ";

    bool first = true;
    for (int i = 0; i <= n; i++) {
        f c = coef[i];
        int power = n - i;
        if (c == 0) continue;

        if (!first) {
            if (c > 0) out << " + ";
            else out << " - ";
        } else {
            if (c < 0) out << "-";
        }

        if (abs(c) != 1 || power == 0)
            out << abs(c);

        if (power > 0) {
            out << "x";
            if (power > 1)
                out << "^" << power;
        }

        first = false;
    }
    out << "\n\n";
}

// Print header
void printHeader(ofstream &out) {
    out << "============================================\n";
    out << "          Newton-Raphson Method\n";
    out << "============================================\n\n";
}
void printHeader() {
    cout << "============================================\n";
    cout << "          Newton-Raphson Method\n";
    cout << "============================================\n\n";
}

// Print roots table
void printRootsTable(ofstream &out, const vector<f> &roots, const vector<f> &coef, f tolerance = 1e-6) {
    out << setw(6) << "Index" 
        << setw(20) << "Root Value" << "\n";
    out << string(26, '-') << "\n";

    out << fixed << setprecision(6);
    for (size_t i = 0; i < roots.size(); i++) {
        f fx = fun(coef, roots[i]);
        out << setw(6) << (i + 1)
            << setw(20) << roots[i]
            << "\n";
    }
    out << "\n";
}

int main() {
    string inputFile, outputFile;

    printHeader();
    // File names
    cout << "Enter input file name: ";
    cin >> inputFile;
    cout << "Enter output file name: ";
    cin >> outputFile;

    ifstream in(inputFile);
    ofstream out(outputFile);

    if (!in.is_open() || !out.is_open()) {
        cout << "Error: Cannot open input/output file!\n";
        return 1;
    }

    int degree;
    in >> degree;
    vector<f> coef(degree + 1);
    for (int i = 0; i <= degree; i++)
        in >> coef[i];

    printHeader(out);

    out << "Polynomial Degree : " << degree << "\n";
    out << "Coefficients      : ";
    for (auto c : coef) out << c << " ";
    out << "\n\n";

    printPolynomial(out, coef);

    // Parameters
    f searchRange = 5000.0, step = 0.5, tolerance = 1e-6;
    vector<f> roots, intervals;

    // Find intervals with sign changes
    for (f i = -searchRange; i <= searchRange; i += step) {
        f fx1 = fun(coef, i);
        f fx2 = fun(coef, i + step);

        // Exact root at grid
        if (abs(fx1) < tolerance) {
            bool duplicate = false;
            for (auto r : roots) if (abs(r - i) < tolerance) { duplicate = true; break; }
            if (!duplicate) roots.push_back(i);
        }
        // Sign change
        else if (fx1 * fx2 < 0.0) {
            intervals.push_back(i);
        }
    }

    // Apply Newton-Raphson Method
    for (f start : intervals) {
        f x = start;
        map<f, bool> visited;

        while (true) {
            f fx = fun(coef, x);
            f dfx = derivative(coef, x);

            // Check for zero derivative (avoid division by zero)
            if (abs(dfx) < 1e-10) break;

            // Newton-Raphson formula: x_new = x - f(x)/f'(x)
            f x_new = x - (fx / dfx);
            f fx_new = fun(coef, x_new);

            // Round for map comparison to avoid floating point precision issues
            f x_new_rounded = round(x_new * 1e6) / 1e6;

            // Check convergence or repeated value
            if (abs(fx_new) < tolerance || visited[x_new_rounded]) {
                bool duplicate = false;
                for (auto r : roots) if (abs(r - x_new) < tolerance) { duplicate = true; break; }
                if (!duplicate) roots.push_back(x_new);
                break;
            }

            visited[x_new_rounded] = true;

            // Update for next iteration
            x = x_new;
        }
    }

    sort(roots.begin(), roots.end());

    // Output roots table
    out << "Roots\n";
    out << "============================================\n";
    if (roots.empty()) {
        out << "No real roots found in the given range.\n\n";
    } else {
        printRootsTable(out, roots, coef, tolerance);
    }

    out << "============================================\n";
    out << "Computation Completed Successfully.\n";

    in.close();
    out.close();

    cout << "Computation completed. Results written to '" << outputFile << "'\n";

    return 0;
}

```

---

## Newton-Raphson Input
**Input1 (input1.txt):**   
```
3
1 -6 11 -6
```
**Input2 (input2.txt):**   
```
2
1 0 1
```
---

## Newton-Raphson Output
**Output1 (output1.txt):** 
```
============================================
          Newton-Raphson Method
============================================

Polynomial Degree : 3
Coefficients      : 1 -6 11 -6 

f(x) = x^3 - 6x^2 + 11x - 6

Roots
============================================
 Index          Root Value
--------------------------
     1            1.000000
     2            2.000000
     3            3.000000

============================================
Computation Completed Successfully.

```

**Output2 (output2.txt):** 
```
============================================
          Newton-Raphson Method
============================================

Polynomial Degree : 2
Coefficients      : 1 0 1 

f(x) = x^2 + 1

Roots
============================================
No real roots found in the given range.

============================================
Computation Completed Successfully.

```
---
# Secant Method
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Solution%20of%20Non-linear%20Equations/Secant%20Method/)

## Secant Theory
### Mathematical Foundation

The Secant Method is based on **linear interpolation** between two points to approximate the root. Unlike bracketing methods, it does not require the initial points to bracket a root, but it benefits from sign changes for better convergence.

**Core Principle:**

Given two points $x_1$ and $x_2$ with corresponding function values $f(x_1)$ and $f(x_2)$, the method approximates the root by finding where the secant line through these points crosses the x-axis.

**Secant Formula:**

$$x_0 = x_1 - f(x_1) \cdot \frac{x_2 - x_1}{f(x_2) - f(x_1)}$$

This can also be written as:

$$x_0 = \frac{x_1 f(x_2) - x_2 f(x_1)}{f(x_2) - f(x_1)}$$

**Geometric Interpretation:**

The Secant Method draws a straight line (secant) connecting the points $(x_1, f(x_1))$ and $(x_2, f(x_2))$, then finds the x-intercept of this line. This x-intercept becomes the new approximation.

**Derivative Approximation:**

The Secant Method approximates the derivative in Newton-Raphson:

$$f'(x) \approx \frac{f(x_2) - f(x_1)}{x_2 - x_1}$$

Substituting this into Newton's formula $x_{\text{new}} = x_{\text{old}} - \frac{f(x_{\text{old}})}{f'(x_{\text{old}})}$ gives the Secant formula.

**Iteration Process:**

1. Start with two initial points: $x_1$ and $x_2$
2. Calculate new approximation using the Secant formula: $x_0 = x_1 - f(x_1) \cdot \frac{x_2 - x_1}{f(x_2) - f(x_1)}$
3. Check convergence: $|f(x_0)| < \epsilon$
4. Update points: $x_1 = x_2$, $x_2 = x_0$ (shift the two most recent points)
5. Repeat until convergence or cycle detection

**Convergence Criteria:**

The algorithm stops when:

$$|f(x_0)| < \epsilon$$

where:
- $\epsilon = 10^{-6}$ (function value tolerance)

Additionally, the implementation includes **cycle detection** using a map to prevent infinite loops when the method enters a repetitive pattern.

**Comparison with Other Methods:**

| Aspect | Bisection | False Position | Secant | Newton-Raphson |
|--------|-----------|----------------|--------|----------------|
| Derivative Required | No | No | No | Yes |
| Convergence Rate | Linear | Superlinear | Superlinear (≈1.618) | Quadratic (2.0) |
| Initial Points | 2 (bracket) | 2 (bracket) | 2 (any) | 1 |
| Bracketing Required | Yes | Yes | No | No |
| Function Evaluations/Iteration | 1 | 2 | 1 | 1 + derivative |
| Convergence Order | 1.0 | ~1.6 (varies) | φ ≈ 1.618 | 2.0 |
| Best For | Guaranteed slow | Nearly linear | General purpose | When derivative available |

**Convergence Order:**

The Secant Method has a convergence order of $\phi = \frac{1 + \sqrt{5}}{2} \approx 1.618$ (the golden ratio), which means:

$$|e_{n+1}| \approx C |e_n|^{1.618}$$

where $e_n$ is the error at iteration $n$. This is faster than linear (order 1) but slower than quadratic (order 2).

### Algorithm Steps

1. **Input Processing**:
   - Read polynomial degree $n$
   - Read $n+1$ coefficients: $a_0, a_1, ..., a_n$ for $f(x) = a_0 x^n + a_1 x^{n-1} + ... + a_{n-1}x + a_n$

2. **Interval Scanning**:
   - Define search range: $[-R, R]$ where $R = 5000$ (default)
   - Use step size: $h = 0.5$
   - For each point $x_i = -R + i \cdot h$:
     - Check if $|f(x_i)| < \epsilon$ (direct root detection)
     - Check if $f(x_i) \cdot f(x_{i+1}) < 0$ (sign change detection)
     - Store intervals with sign changes

3. **Secant Iteration** (for each interval):
   - Initialize: $x_1 = x_i$, $x_2 = x_{i+1}$ (two consecutive points)
   - Create empty visited map for cycle detection
   - While not converged:
     1. Evaluate: $f(x_1) = f(x_1)$, $f(x_2) = f(x_2)$
     2. **Check for division by zero:** If $|f(x_2) - f(x_1)| < 10^{-10}$, break
     3. **Apply Secant formula:** 
        $$x_0 = x_1 - f(x_1) \cdot \frac{x_2 - x_1}{f(x_2) - f(x_1)}$$
     4. Evaluate: $f_0 = f(x_0)$
     5. Round $x_0$ for map: $x_{0,\text{rounded}} = \text{round}(x_0 \times 10^6) / 10^6$
     6. **Convergence/cycle check:**
        - If $|f_0| < \epsilon$ OR $x_{0,\text{rounded}}$ is in visited map:
          - Check for duplicate in roots list
          - Add $x_0$ to roots if unique
          - Break loop
     7. Mark $x_{0,\text{rounded}}$ as visited
     8. **Update for next iteration:** $x_1 = x_2$, $x_2 = x_0$ (shift points)

4. **Root Collection**:
   - Combine roots from direct detection and Secant iterations
   - Filter duplicates (within tolerance $\epsilon$)
   - Sort roots in ascending order

5. **Output Generation**:
   - Display polynomial equation
   - Print final roots table with indices
   - Report completion status

### Complexity Analysis

- **Time Complexity**:
  - Interval scanning: $O(R/h)$ where $R$ is search range, $h$ is step size
  - Single root finding: $O(k)$ where $k$ is number of iterations
  - For $m$ roots: $O(R/h + m \cdot k)$
  - Typically requires **40-60% fewer iterations** than Bisection

- **Space Complexity**: 
  - $O(m + k)$ where $m$ is number of roots
  - $k$ is average size of visited map per root (cycle detection overhead)

- **Convergence Rate**: 
  - **Superlinear** with order $\phi \approx 1.618$
  - Error reduces by factor of approximately $e^{1.618}$ each iteration
  - Significantly faster than Bisection (linear) and False Position
  - Slightly slower than Newton-Raphson (quadratic) but no derivative needed

- **Function Evaluations**:
  - **2 evaluations per iteration** ($f(x_1)$ and $f(x_2)$)
  - More efficient than Newton-Raphson when derivative computation is expensive

- **Accuracy**: 
  - Guaranteed to converge to a root within tolerance $\epsilon$ (when convergent)
  - Final error: $|r_{\text{true}} - r_{\text{computed}}| \lesssim \epsilon$
  - May not converge if initial points are poorly chosen (far from root)
---

## Secant Code
```cpp
#include <bits/stdc++.h>
using namespace std;

#define f double

// Polynomial function f(x)
f fun(const vector<f> &coef, f x) {
    f result = 0.0;
    int n = coef.size() - 1;
    for (int i = 0; i <= n; i++)
        result += coef[i] * pow(x, n - i);
    return result;
}

// Print polynomial in readable form
void printPolynomial(ofstream &out, const vector<f> &coef) {
    int n = coef.size() - 1;
    out << "f(x) = ";

    bool first = true;
    for (int i = 0; i <= n; i++) {
        f c = coef[i];
        int power = n - i;
        if (c == 0) continue;

        if (!first) {
            if (c > 0) out << " + ";
            else out << " - ";
        } else {
            if (c < 0) out << "-";
        }

        if (abs(c) != 1 || power == 0)
            out << abs(c);

        if (power > 0) {
            out << "x";
            if (power > 1)
                out << "^" << power;
        }

        first = false;
    }
    out << "\n\n";
}

// Print header
void printHeader(ofstream &out) {
    out << "============================================\n";
    out << "              Secant Method\n";
    out << "============================================\n\n";
}
void printHeader() {
    cout << "============================================\n";
    cout << "              Secant Method\n";
    cout << "============================================\n\n";
}

// Print roots table
void printRootsTable(ofstream &out, const vector<f> &roots, const vector<f> &coef, f tolerance = 1e-6) {
    out << setw(6) << "Index" 
        << setw(20) << "Root Value" << "\n";
    out << string(26, '-') << "\n";

    out << fixed << setprecision(6);
    for (size_t i = 0; i < roots.size(); i++) {
        f fx = fun(coef, roots[i]);
        out << setw(6) << (i + 1)
            << setw(20) << roots[i]
            << "\n";
    }
    out << "\n";
}

int main() {
    string inputFile, outputFile;

    printHeader();
    // File names
    cout << "Enter input file name: ";
    cin >> inputFile;
    cout << "Enter output file name: ";
    cin >> outputFile;

    ifstream in(inputFile);
    ofstream out(outputFile);

    if (!in.is_open() || !out.is_open()) {
        cout << "Error: Cannot open input/output file!\n";
        return 1;
    }

    int degree;
    in >> degree;
    vector<f> coef(degree + 1);
    for (int i = 0; i <= degree; i++)
        in >> coef[i];

    printHeader(out);

    out << "Polynomial Degree : " << degree << "\n";
    out << "Coefficients      : ";
    for (auto c : coef) out << c << " ";
    out << "\n\n";

    printPolynomial(out, coef);

    // Parameters
    f searchRange = 5000.0, step = 0.5, tolerance = 1e-6;
    vector<f> roots, intervals;

    // Find intervals with sign changes
    for (f i = -searchRange; i <= searchRange; i += step) {
        f fx1 = fun(coef, i);
        f fx2 = fun(coef, i + step);

        // Exact root at grid
        if (abs(fx1) < tolerance) {
            bool duplicate = false;
            for (auto r : roots) if (abs(r - i) < tolerance) { duplicate = true; break; }
            if (!duplicate) roots.push_back(i);
        }
        // Sign change
        else if (fx1 * fx2 < 0.0) {
            intervals.push_back(i);
        }
    }

    // Apply Secant Method
    for (f start : intervals) {
        f x1 = start, x2 = start + step;
        map<f, bool> visited;

        while (true) {
            f fx1 = fun(coef, x1);
            f fx2 = fun(coef, x2);

            // Secant formula
            if (abs(fx2 - fx1) < 1e-10) break; // Avoid division by zero

            f x0 = x1 - fx1 * ((x2 - x1) / (fx2 - fx1));
            f fx0 = fun(coef, x0);

            // Round for map comparison to avoid floating point precision issues
            f x0_rounded = round(x0 * 1e6) / 1e6;

            // Check convergence or repeated value
            if (abs(fx0) < tolerance || visited[x0_rounded]) {
                bool duplicate = false;
                for (auto r : roots) if (abs(r - x0) < tolerance) { duplicate = true; break; }
                if (!duplicate) roots.push_back(x0);
                break;
            }

            visited[x0_rounded] = true;

            // Update for next iteration (Secant method updates)
            x1 = x2;
            x2 = x0;
        }
    }

    sort(roots.begin(), roots.end());

    // Output roots table
    out << "Roots\n";
    out << "============================================\n";
    if (roots.empty()) {
        out << "No real roots found in the given range.\n\n";
    } else {
        printRootsTable(out, roots, coef, tolerance);
    }

    out << "============================================\n";
    out << "Computation Completed Successfully.\n";

    in.close();
    out.close();

    cout << "Computation completed. Results written to '" << outputFile << "'\n";

    return 0;
}

```

---

## Secant Input
**Input1 (input1.txt):**   
```
3
1 -6 11 -6
```
**Input2 (input2.txt):**   
```
2
1 0 1

```
---

## Secant Output
**Output1 (output1.txt):** 
```
============================================
              Secant Method
============================================

Polynomial Degree : 3
Coefficients      : 1 -6 11 -6 

f(x) = x^3 - 6x^2 + 11x - 6

Roots
============================================
 Index          Root Value
--------------------------
     1            1.000000
     2            2.000000
     3            3.000000

============================================
Computation Completed Successfully.

```

**Output2 (output2.txt):** 
```
============================================
              Secant Method
============================================

Polynomial Degree : 2
Coefficients      : 1 0 1 

f(x) = x^2 + 1

Roots
============================================
No real roots found in the given range.

============================================
Computation Completed Successfully.

```
---
---
# Interpolation and Approximation
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-cyan?style=for-the-badge)](./Interpolation%20and%20Approximation/)
## 📐 Mathematical Foundation
### What is Interpolation?

**Interpolation** is the process of constructing new data points within the range of a discrete set of known data points. 

**Formal Definition**: Given n+1 distinct points `(x₀, y₀), (x₁, y₁), ..., (xₙ, yₙ)`, find a function f(x) from a specified class of functions such that: 

```
f(xᵢ) = yᵢ  for i = 0, 1, 2, ..., n
```

**Key Properties**:
- **Interpolation constraint**: The function must pass through all given data points
- **Uniqueness**: For polynomial interpolation of degree ≤ n, the solution is unique
- **Continuity**: The interpolating function is continuous (unlike step functions)

**Interpolation vs.  Extrapolation**:
```
Data points:   x₀  x₁  x₂  x₃  x₄
               |---|---|---|---|
              
Interpolation:  Estimate within [x₀, x₄]  ✓ More reliable
Extrapolation: Estimate outside [x₀, x₄] ⚠️  Less reliable, higher error
```

### What is Approximation?

**Approximation** is the broader process of finding a function that is "close to" the data points, but doesn't necessarily pass through them exactly.

**Key Differences**: 

| Aspect | Interpolation | Approximation |
|--------|--------------|---------------|
| **Exactness** | f(xᵢ) = yᵢ (exact) | f(xᵢ) ≈ yᵢ (approximate) |
| **Use Case** | Exact data | Noisy data, trend fitting |
| **Flexibility** | Must pass through points | Can smooth out noise |
| **Degree** | Degree n for n+1 points | Usually lower degree |
| **Examples** | Newton, Lagrange | Least squares, splines |

**When to Use Each**:
- **Interpolation**: Data is accurate, need exact fit
- **Approximation**: Data has noise/errors, want smooth trend

### Polynomial Interpolation

**The Fundamental Theorem**:  Given n+1 distinct points, there exists a **unique polynomial of degree ≤ n** that passes through all points.

**General Form**: 
```
Pₙ(x) = a₀ + a₁x + a₂x² + ...  + aₙxⁿ
```

**Why Polynomials?**  
✅ Easy to compute and differentiate  
✅ Continuous and smooth  
✅ Well-understood mathematical properties  
✅ Can approximate many functions (Weierstrass Approximation Theorem)  

**Challenges**:  
❌ **Runge's Phenomenon**: High-degree polynomials can oscillate wildly between points  
❌ **Sensitivity**: Small changes in data can cause large changes in polynomial  
❌ **Computational Cost**: Higher degrees require more computation  

### Fundamental Concepts

#### **Finite Differences**

Finite differences measure how function values change between consecutive points.  They're the discrete analog of derivatives.

**Forward Difference Operator (Δ)**:
```
Δy₀ = y₁ - y₀
Δ²y₀ = Δy₁ - Δy₀ = (y₂ - y₁) - (y₁ - y₀) = y₂ - 2y₁ + y₀
Δⁿyᵢ = Δⁿ⁻¹yᵢ₊₁ - Δⁿ⁻¹yᵢ
```

**Backward Difference Operator (∇)**:
```
∇yₙ = yₙ - yₙ₋₁
∇²yₙ = ∇yₙ - ∇yₙ₋₁ = (yₙ - yₙ₋₁) - (yₙ₋₁ - yₙ₋₂) = yₙ - 2yₙ₋₁ + yₙ₋₂
∇ⁿyᵢ = ∇ⁿ⁻¹yᵢ - ∇ⁿ⁻¹yᵢ₋₁
```

**Divided Difference** (for non-uniform spacing):
```
f[x₀, x₁] = (f(x₁) - f(x₀)) / (x₁ - x₀)
f[x₀, x₁, x₂] = (f[x₁, x₂] - f[x₀, x₁]) / (x₂ - x₀)
```

**Visual Example**: Forward Difference Table for y = x²

```
x    y    Δy   Δ²y   Δ³y
0    0    
           1
1    1          2
           3         0
2    4          2
           5         0
3    9          2
           7
4    16
```

#### **Difference Tables**

A **difference table** organizes successive differences in a triangular structure, making patterns visible and computation efficient.

**Properties**:
- For a polynomial of degree n, the nth differences are constant
- Higher-order differences become zero
- Helps detect errors in data (irregular patterns indicate problems)

#### **Interpolating Polynomial Uniqueness**

**Theorem**: Given n+1 distinct points, the interpolating polynomial of degree ≤ n is unique.

**Proof Sketch**:
- Suppose P(x) and Q(x) both interpolate the data
- Then R(x) = P(x) - Q(x) has n+1 roots (at all data points)
- But R(x) has degree ≤ n
- A polynomial of degree n can have at most n roots (unless it's zero)
- Therefore R(x) ≡ 0, so P(x) = Q(x)

**Implication**: All interpolation methods (Newton, Lagrange, etc.) produce the same polynomial, just in different forms! 

#### **Error Analysis**

The **interpolation error** measures how far the interpolating polynomial is from the true function.

**Error Bound** (for polynomial interpolation):
```
|f(x) - Pₙ(x)| ≤ (Mₙ₊₁/(n+1)!) × |(x-x₀)(x-x₁)...(x-xₙ)|
```

Where Mₙ₊₁ = max|f⁽ⁿ⁺¹⁾(x)| over the interval. 

**Key Insights**:
- Error depends on f⁽ⁿ⁺¹⁾ (smoothness of function)
- Error is smallest near the middle of data points
- Error grows near the boundaries (especially for high-degree polynomials)
- Product term |(x-x₀)(x-x₁)...(x-xₙ)| is zero at data points (as expected)
---
# Newton's Forward Interpolation
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Interpolation%20and%20Approximation/Newton's%20Forward%20Interpolation/)

## Forward Interpolation Theory
### Mathematical Foundation

Given a set of data points $(x_0, y_0), (x_1, y_1), ..., (x_n, y_n)$ that are equally spaced with step size $h = x_{i+1} - x_i$, Newton's Forward Interpolation formula is:

$$P(x) = y_0 + u \Delta y_0 + \frac{u(u-1)}{2!} \Delta^2 y_0 + \frac{u(u-1)(u-2)}{3!} \Delta^3 y_0 + ...  + \frac{u(u-1).. .(u-n+1)}{n!} \Delta^n y_0$$

where:  
- $u = \frac{x - x_0}{h}$ (normalized position)
- $\Delta^k y_0$ is the $k$-th forward difference at $x_0$

**Forward Difference Operator:**

The forward differences are computed recursively:
- $\Delta y_i = y_{i+1} - y_i$ (first forward difference)
- $\Delta^2 y_i = \Delta y_{i+1} - \Delta y_i$ (second forward difference)
- $\Delta^k y_i = \Delta^{k-1} y_{i+1} - \Delta^{k-1} y_i$ (k-th forward difference)

**Forward Difference Table:**

| $x$ | $y$ | $\Delta y$ | $\Delta^2 y$ | $\Delta^3 y$ | $\Delta^4 y$ |
|-----|-----|-----------|-------------|-------------|-------------|
| $x_0$ | $y_0$ | $\Delta y_0$ | $\Delta^2 y_0$ | $\Delta^3 y_0$ | $\Delta^4 y_0$ |
| $x_1$ | $y_1$ | $\Delta y_1$ | $\Delta^2 y_1$ | $\Delta^3 y_1$ | |
| $x_2$ | $y_2$ | $\Delta y_2$ | $\Delta^2 y_2$ | | |
| $x_3$ | $y_3$ | $\Delta y_3$ | | | |
| $x_4$ | $y_4$ | | | | |

### Algorithm Steps

1. **Verify Equal Spacing**:  Check that all data points have uniform spacing $h$
2. **Build Forward Difference Table**: 
   - Initialize first column with $y$ values
   - Compute successive forward differences using:  $\Delta^k y_i = \Delta^{k-1} y_{i+1} - \Delta^{k-1} y_i$
3. **Calculate Normalized Position**: $u = \frac{x - x_0}{h}$
4. **Apply Newton's Forward Formula**:
   - Start with $P(x) = y_0$
   - Add terms:  $\frac{u(u-1)...(u-k+1)}{k!} \Delta^k y_0$ for $k = 1, 2, ..., n$
5. **Return Interpolated Value**

### Complexity Analysis

- **Time Complexity**:
  - Building difference table: $O(n^2)$ where $n$ is the number of data points
  - Single interpolation: $O(n)$
  - Total for $m$ interpolations: $O(n^2 + mn)$
  
- **Space Complexity**: $O(n^2)$ for storing the forward difference table

- **Accuracy**: For equally spaced points, Newton's Forward Interpolation provides exact results for polynomial data of degree $\leq n-1$
---

## Forward Interpolation Code
```cpp
#include <bits/stdc++.h>
using namespace std;

/*
   Utility:  factorial up to n (double)
*/
double factorial(int n){
    double f = 1.0;
    for (int i=2; i<=n; i++) f *= i;
    return f;
}

/*
   Print Data Points Table
*/
void printDataTable(const vector<double>& xs, const vector<double>& ys, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  DATA POINTS TABLE\n";
    out << "====================================\n";
    
    // Print x values
    out << setw(10) << "x";
    for (int i=0; i<n; i++){
        out << setw(15) << fixed << setprecision(4) << xs[i];
    }
    out << "\n";
    
    // Print separator
    out << string(10 + 15*n, '-') << "\n";
    
    // Print y values
    out << setw(10) << "y";
    for (int i=0; i<n; i++){
        out << setw(15) << fixed << setprecision(6) << ys[i];
    }
    out << "\n";
    out << "====================================\n";
}

/*
   Build Forward Difference Table
*/
vector<vector<double>> buildForwardDiffTable(const vector<double>& xs, const vector<double>& ys) {
    int n = (int)xs.size();
    vector<vector<double>> diff(n, vector<double>(n));
    
    for (int i=0; i<n; i++) diff[i][0] = ys[i];
    for (int j=1; j<n; j++){
        for (int i=0; i<n-j; i++){
            diff[i][j] = diff[i+1][j-1] - diff[i][j-1];
        }
    }
    return diff;
}

/*
   Print Forward Difference Table
*/
void printForwardDiffTable(const vector<double>& xs, const vector<vector<double>>& diff, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  FORWARD DIFFERENCE TABLE\n";
    out << "====================================\n";
    out << setw(12) << "x" << setw(15) << "y";
    for (int j=1; j<n; j++){
        out << setw(15) << ("Delta^" + to_string(j) + "y");
    }
    out << "\n" << string(12 + 15*n, '-') << "\n";
    
    for (int i=0; i<n; i++){
        out << fixed << setprecision(4) << setw(12) << xs[i];
        for (int j=0; j<n-i; j++){
            out << setprecision(6) << setw(15) << diff[i][j];
        }
        out << "\n";
    }
    out << "====================================\n";
}

/*
   Newton's Forward Interpolation using pre-built table
*/
double newtonForwardWithTable(const vector<double>& xs, const vector<vector<double>>& diff, double x) {
    int n = (int)xs.size();
    if (n == 1) return diff[0][0];
    
    double h = xs[1] - xs[0];
    double u = (x - xs[0]) / h;
    double result = diff[0][0];
    double u_prod = 1.0;
    
    for (int k=1; k<n; k++){
        u_prod *= (u - (k-1));
        result += (u_prod / factorial(k)) * diff[0][k];
    }
    return result;
}

/*
   Process and output interpolation results
*/
void processInterpolation(const vector<double>& xs, const vector<vector<double>>& diff,
                         const vector<double>& xInterpolate, vector<double>& results,
                         ostream& cout_stream, ostream& fout) {
    int m = xInterpolate.size();
    int n = xs.size();
    
    for (int i=0; i<m; i++){
        double result = newtonForwardWithTable(xs, diff, xInterpolate[i]);
        results[i] = result;
        
        bool isExtrap = (xInterpolate[i] < xs[0] || xInterpolate[i] > xs[n-1]);
        string extrapNote = isExtrap ?  " (Extrapolation)" : "";
        
        string output = "Point " + to_string(i+1) + ": x = ";
        cout_stream << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        cout_stream << "         y = " << setprecision(6) << result << "\n";
        
        fout << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        fout << "         y = " << setprecision(6) << result << "\n";
    }
}

/*
   Print section header
*/
void printHeader(const string& title, ostream& out) {
    out << "\n====================================\n";
    out << title << "\n";
    out << "====================================\n";
}

/*
   Main with File I/O
*/
int main() {
    printHeader("NEWTON'S FORWARD INTERPOLATION", cout);
    
    string inputFile, outputFile;
    cout << "\nEnter input filename: ";
    cin >> inputFile;
    cout << "Enter output filename: ";
    cin >> outputFile;
    
    // Open input file
    ifstream fin(inputFile);
    if (!fin) {
        cerr << "Error: Cannot open input file '" << inputFile << "'\n";
        return 1;
    }
    
    // Read data points
    int n;
    fin >> n;
    if (n <= 0) {
        cerr << "Error: Invalid number of data points\n";
        return 1;
    }
    
    vector<double> xs(n), ys(n);
    for (int i=0; i<n; i++){
        fin >> xs[i] >> ys[i];
    }
    
    // Read interpolation points
    int m;
    fin >> m;
    if (m <= 0) {
        cerr << "Error: Invalid number of interpolation points\n";
        return 1;
    }
    
    vector<double> xInterpolate(m);
    for (int i=0; i<m; i++){
        fin >> xInterpolate[i];
    }
    
    // Check for additional data point
    double xAdditional, yAdditional;
    bool hasAdditional = false;
    if (fin >> xAdditional >> yAdditional) {
        hasAdditional = true;
    }
    fin.close();
    
    // Verify equal spacing
    if (n > 1) {
        double h = xs[1] - xs[0];
        for (int i=2; i<n; i++){
            if (fabs((xs[i] - xs[i-1]) - h) > 1e-9) {
                cerr << "Error:  Data points are not equally spaced!\n";
                return 1;
            }
        }
        cout << "\nNumber of data points: " << n << "\n";
        cout << "Step size (h): " << fixed << setprecision(6) << h << "\n";
    }
    
    // Open output file
    ofstream fout(outputFile);
    if (!fout) {
        cerr << "Error:  Cannot create output file '" << outputFile << "'\n";
        return 1;
    }
    
    printHeader("NEWTON'S FORWARD INTERPOLATION", fout);
    fout << "\nNumber of data points: " << n << "\n";
    if (n > 1) {
        fout << "Step size (h): " << fixed << setprecision(6) << (xs[1] - xs[0]) << "\n";
    }
    
    // Print data points table
    printDataTable(xs, ys, cout);
    printDataTable(xs, ys, fout);
    
    // Build and print forward difference table
    vector<vector<double>> diffTable = buildForwardDiffTable(xs, ys);
    printForwardDiffTable(xs, diffTable, cout);
    printForwardDiffTable(xs, diffTable, fout);
    
    // Perform interpolation
    printHeader("  INTERPOLATION RESULTS", cout);
    printHeader("  INTERPOLATION RESULTS", fout);
    
    vector<double> interpolatedValues(m);
    processInterpolation(xs, diffTable, xInterpolate, interpolatedValues, cout, fout);
    
    cout << "====================================\n";
    fout << "====================================\n";
    
    // Process additional point if provided
    if (hasAdditional) {
        vector<double> xsNew = xs;
        vector<double> ysNew = ys;
        
        // Insert in sorted order
        int insertPos = n;
        for (int i=0; i<n; i++){
            if (xAdditional < xs[i]) {
                insertPos = i;
                break;
            }
        }
        xsNew.insert(xsNew.begin() + insertPos, xAdditional);
        ysNew.insert(ysNew. begin() + insertPos, yAdditional);
        
        // Verify equal spacing
        int nNew = n + 1;
        double hNew = xsNew[1] - xsNew[0];
        bool validSpacing = true;
        for (int i=2; i<nNew; i++){
            if (fabs((xsNew[i] - xsNew[i-1]) - hNew) > 1e-9) {
                validSpacing = false;
                break;
            }
        }
        
        if (validSpacing) {
            printHeader("  WITH ADDITIONAL DATA POINT", cout);
            printHeader("  WITH ADDITIONAL DATA POINT", fout);
            
            cout << "Additional point:  x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            cout << "New number of data points: " << nNew << "\n";
            cout << "New step size (h): " << setprecision(6) << hNew << "\n";
            
            fout << "Additional point: x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            fout << "New number of data points: " << nNew << "\n";
            fout << "New step size (h): " << setprecision(6) << hNew << "\n";
            
            // Print updated data points table
            printDataTable(xsNew, ysNew, cout);
            printDataTable(xsNew, ysNew, fout);
            
            // Build new difference table
            vector<vector<double>> diffTableNew = buildForwardDiffTable(xsNew, ysNew);
            printForwardDiffTable(xsNew, diffTableNew, cout);
            printForwardDiffTable(xsNew, diffTableNew, fout);
            
            // Re-interpolate and show differences
            printHeader("  UPDATED INTERPOLATION RESULTS", cout);
            printHeader("  UPDATED INTERPOLATION RESULTS", fout);
            
            for (int i=0; i<m; i++){
                double resultOld = interpolatedValues[i];
                double resultNew = newtonForwardWithTable(xsNew, diffTableNew, xInterpolate[i]);
                double absError = fabs(resultNew - resultOld);
                double relError = (fabs(resultNew) > 1e-15) ? (absError / fabs(resultNew)) * 100.0 : 0.0;
                
                cout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                cout << "         Old y = " << setprecision(6) << resultOld << "\n";
                cout << "         New y = " << setprecision(6) << resultNew << "\n";
                cout << "         Absolute Difference:   " << scientific << setprecision(6) << absError << "\n";
                cout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
                
                fout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                fout << "         Old y = " << setprecision(6) << resultOld << "\n";
                fout << "         New y = " << setprecision(6) << resultNew << "\n";
                fout << "         Absolute Difference:  " << scientific << setprecision(6) << absError << "\n";
                fout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
            }
            
            cout << "====================================\n";
            fout << "====================================\n";
        } else {
            printHeader("WARNING: Additional point breaks equal spacing!", cout);
            cout << "Cannot use Newton's Forward Interpolation.\n";
            cout << "====================================\n";
            
            printHeader("WARNING: Additional point breaks equal spacing!", fout);
            fout << "Cannot use Newton's Forward Interpolation.\n";
            fout << "====================================\n";
        }
    }
    
    fout.close();
    cout << "\nResults have been written to '" << outputFile << "'\n";
    cout << "Program executed successfully!\n";
    
    return 0;
}
```

---

## Forward Interpolation Input
### Input File Format

The input file must follow this specific structure: 

```
n                    # Line 1: Number of data points (integer)
x₁ y₁               # Line 2: First data point (two space-separated numbers)
x₂ y₂               # Line 3: Second data point
x₃ y₃               # Line 4: Third data point
...
xₙ yₙ               # Line n+1: n-th data point
m                   # Line n+2: Number of interpolation points (integer)
x_interp₁           # Line n+3: First x-value to interpolate (single number)
x_interp₂           # Line n+4: Second x-value to interpolate
...
x_interpₘ           # Line n+2+m: m-th x-value to interpolate
x_add y_add         # Line n+3+m (OPTIONAL): Additional data point for error analysis
```
**Input1 (input1.txt):**   
```
5
0.0 1.0
0.5 1.6487
1.0 2.7183
1.5 4.4817
2.0 7.3891
3
0.25
0.75
1.25
```
**Input2 (input2.txt):**   
```
5
0.0 1.0
0.5 1.6487
1.0 2.7183
1.5 4.4817
2.0 7.3891
3
0.25
0.75
1.25
2.5 12.1825
```
---

## Forward Interpolation Output
**Output1 (output1.txt):** 
```

====================================
NEWTON'S FORWARD INTERPOLATION
====================================

Number of data points: 5
Step size (h): 0.500000

====================================
  DATA POINTS TABLE
====================================
         x         0.0000         0.5000         1.0000         1.5000         2.0000
-------------------------------------------------------------------------------------
         y       1.000000       1.648700       2.718300       4.481700       7.389100
====================================

====================================
  FORWARD DIFFERENCE TABLE
====================================
           x              y       Delta^1y       Delta^2y       Delta^3y       Delta^4y
---------------------------------------------------------------------------------------
      0.0000       1.000000       0.648700       0.420900       0.272900       0.177300
      0.5000       1.648700       1.069600       0.693800       0.450200
      1.0000       2.718300       1.763400       1.144000
      1.5000       4.481700       2.907400
      2.0000       7.389100
====================================

====================================
  INTERPOLATION RESULTS
====================================
Point 1: x = 0.250000
         y = 1.281868
Point 2: x = 0.750000
         y = 2.117987
Point 3: x = 1.250000
         y = 3.489293
====================================

```

**Output2 (output2.txt):** 
```

====================================
NEWTON'S FORWARD INTERPOLATION
====================================

Number of data points: 5
Step size (h): 0.500000

====================================
  DATA POINTS TABLE
====================================
         x         0.0000         0.5000         1.0000         1.5000         2.0000
-------------------------------------------------------------------------------------
         y       1.000000       1.648700       2.718300       4.481700       7.389100
====================================

====================================
  FORWARD DIFFERENCE TABLE
====================================
           x              y       Delta^1y       Delta^2y       Delta^3y       Delta^4y
---------------------------------------------------------------------------------------
      0.0000       1.000000       0.648700       0.420900       0.272900       0.177300
      0.5000       1.648700       1.069600       0.693800       0.450200
      1.0000       2.718300       1.763400       1.144000
      1.5000       4.481700       2.907400
      2.0000       7.389100
====================================

====================================
  INTERPOLATION RESULTS
====================================
Point 1: x = 0.250000
         y = 1.281868
Point 2: x = 0.750000
         y = 2.117987
Point 3: x = 1.250000
         y = 3.489293
====================================

====================================
  WITH ADDITIONAL DATA POINT
====================================
Additional point: x = 2.500000, y = 12.182500
New number of data points: 6
New step size (h): 0.500000

====================================
  DATA POINTS TABLE
====================================
         x         0.0000         0.5000         1.0000         1.5000         2.0000         2.5000
----------------------------------------------------------------------------------------------------
         y       1.000000       1.648700       2.718300       4.481700       7.389100      12.182500
====================================

====================================
  FORWARD DIFFERENCE TABLE
====================================
           x              y       Delta^1y       Delta^2y       Delta^3y       Delta^4y       Delta^5y
------------------------------------------------------------------------------------------------------
      0.0000       1.000000       0.648700       0.420900       0.272900       0.177300       0.114500
      0.5000       1.648700       1.069600       0.693800       0.450200       0.291800
      1.0000       2.718300       1.763400       1.144000       0.742000
      1.5000       4.481700       2.907400       1.886000
      2.0000       7.389100       4.793400
      2.5000      12.182500
====================================

====================================
  UPDATED INTERPOLATION RESULTS
====================================
Point 1: x = 0.250000
         Old y = 1.281868
         New y = 1.284999
         Absolute Difference:  3.130859e-03
         Relative Difference: 0.2436%

Point 2: x = 0.750000
         Old y = 2.117987
         New y = 2.116645
         Absolute Difference:  1.341797e-03
         Relative Difference: 0.0634%

Point 3: x = 1.250000
         Old y = 3.489293
         New y = 3.490635
         Absolute Difference:  1.341797e-03
         Relative Difference: 0.0384%

====================================

```
---

# Newton's Backward Interpolation
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Interpolation%20and%20Approximation/Newton's%20Backward%20Interpolation/)

## Backward Interpolation Theory
### Mathematical Foundation

Given a set of data points $(x_0, y_0), (x_1, y_1), ..., (x_n, y_n)$ that are equally spaced with step size $h = x_{i+1} - x_i$, Newton's Backward Interpolation formula is:

$$P(x) = y_n + v \nabla y_n + \frac{v(v+1)}{2!} \nabla^2 y_n + \frac{v(v+1)(v+2)}{3!} \nabla^3 y_n + ...  + \frac{v(v+1).. .(v+n-1)}{n!} \nabla^n y_n$$

where: 
- $v = \frac{x - x_n}{h}$ (normalized position from the end)
- $\nabla^k y_n$ is the $k$-th backward difference at $x_n$

**Backward Difference Operator:**

The backward differences are computed recursively:
- $\nabla y_i = y_i - y_{i-1}$ (first backward difference)
- $\nabla^2 y_i = \nabla y_i - \nabla y_{i-1}$ (second backward difference)
- $\nabla^k y_i = \nabla^{k-1} y_i - \nabla^{k-1} y_{i-1}$ (k-th backward difference)

**Backward Difference Table:**

| $x$ | $y$ | $\nabla y$ | $\nabla^2 y$ | $\nabla^3 y$ | $\nabla^4 y$ |
|-----|-----|-----------|-------------|-------------|-------------|
| $x_0$ | $y_0$ | | | | |
| $x_1$ | $y_1$ | $\nabla y_1$ | | | |
| $x_2$ | $y_2$ | $\nabla y_2$ | $\nabla^2 y_2$ | | |
| $x_3$ | $y_3$ | $\nabla y_3$ | $\nabla^2 y_3$ | $\nabla^3 y_3$ | |
| $x_4$ | $y_4$ | $\nabla y_4$ | $\nabla^2 y_4$ | $\nabla^3 y_4$ | $\nabla^4 y_4$ |

### Algorithm Steps

1. **Verify Equal Spacing**: Check that all data points have uniform spacing $h$
2. **Build Backward Difference Table**:
   - Initialize first column with $y$ values
   - Compute successive backward differences using:  $\nabla^k y_i = \nabla^{k-1} y_i - \nabla^{k-1} y_{i-1}$
3. **Calculate Normalized Position**: $v = \frac{x - x_n}{h}$
4. **Apply Newton's Backward Formula**:
   - Start with $P(x) = y_n$
   - Add terms: $\frac{v(v+1).. .(v+k-1)}{k!} \nabla^k y_n$ for $k = 1, 2, ..., n$
5. **Return Interpolated Value**

### Complexity Analysis

- **Time Complexity**:
  - Building difference table: $O(n^2)$ where $n$ is the number of data points
  - Single interpolation: $O(n)$
  - Total for $m$ interpolations: $O(n^2 + mn)$
  
- **Space Complexity**:  $O(n^2)$ for storing the backward difference table

- **Accuracy**:  For equally spaced points, Newton's Backward Interpolation provides exact results for polynomial data of degree $\leq n-1$, with best accuracy near the end of the data range
---

## Backward Interpolation Code
```cpp
#include <bits/stdc++.h>
using namespace std;

/*
   Utility:   factorial up to n (double)
*/
double factorial(int n){
    double f = 1.0;
    for (int i=2; i<=n; i++) f *= i;
    return f;
}

/*
   Print Data Points Table
*/
void printDataTable(const vector<double>& xs, const vector<double>& ys, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  DATA POINTS TABLE\n";
    out << "====================================\n";
    
    // Print x values
    out << setw(5) << "x" << "\t|";
    for (int i=0; i<n; i++){
        out << setw(10) << fixed << setprecision(4) << xs[i] << "\t|";
    }
    out << "\n";
    
    // Print separator
    out << string(5 + 13*n, '-') << "\n";
    
    // Print y values
    out << setw(5) << "y" << "\t|";
    for (int i=0; i<n; i++){
        out << setw(10) << fixed << setprecision(6) << ys[i] << "\t|";
    }
    out << "\n";
    out << "====================================\n";
}

/*
   Build Backward Difference Table
*/
vector<vector<double>> buildBackwardDiffTable(const vector<double>& xs, const vector<double>& ys) {
    int n = (int)xs.size();
    vector<vector<double>> diff(n, vector<double>(n));
    
    for (int i=0; i<n; i++) diff[i][0] = ys[i];
    for (int j=1; j<n; j++){
        for (int i=n-1; i>=j; i--){
            diff[i][j] = diff[i][j-1] - diff[i-1][j-1];
        }
    }
    return diff;
}

/*
   Print Backward Difference Table
*/
void printBackwardDiffTable(const vector<double>& xs, const vector<vector<double>>& diff, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  BACKWARD DIFFERENCE TABLE\n";
    out << "====================================\n";
    out << setw(12) << "x" << setw(15) << "y";
    for (int j=1; j<n; j++){
        out << setw(15) << ("Nabla^" + to_string(j) + "y");
    }
    out << "\n" << string(12 + 15*n, '-') << "\n";
    
    for (int i=0; i<n; i++){
        out << fixed << setprecision(4) << setw(12) << xs[i];
        for (int j=0; j<=i && j<n; j++){
            out << setprecision(6) << setw(15) << diff[i][j];
        }
        out << "\n";
    }
    out << "====================================\n";
}

/*
   Newton's Backward Interpolation using pre-built table
*/
double newtonBackwardWithTable(const vector<double>& xs, const vector<vector<double>>& diff, double x) {
    int n = (int)xs.size();
    if (n == 1) return diff[0][0];
    
    double h = xs[1] - xs[0];
    double v = (x - xs[n-1]) / h;
    double result = diff[n-1][0];
    double v_prod = 1.0;
    
    for (int k=1; k<n; k++){
        v_prod *= (v + (k-1));
        result += (v_prod / factorial(k)) * diff[n-1][k];
    }
    return result;
}

/*
   Process and output interpolation results
*/
void processInterpolation(const vector<double>& xs, const vector<vector<double>>& diff,
                         const vector<double>& xInterpolate, vector<double>& results,
                         ostream& cout_stream, ostream& fout) {
    int m = xInterpolate.size();
    int n = xs.size();
    
    for (int i=0; i<m; i++){
        double result = newtonBackwardWithTable(xs, diff, xInterpolate[i]);
        results[i] = result;
        
        bool isExtrap = (xInterpolate[i] < xs[0] || xInterpolate[i] > xs[n-1]);
        string extrapNote = isExtrap ? " (Extrapolation)" : "";
        
        string output = "Point " + to_string(i+1) + ": x = ";
        cout_stream << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        cout_stream << "         y = " << setprecision(6) << result << "\n";
        
        fout << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        fout << "         y = " << setprecision(6) << result << "\n";
    }
}

/*
   Print section header
*/
void printHeader(const string& title, ostream& out) {
    out << "\n====================================\n";
    out << title << "\n";
    out << "====================================\n";
}

/*
   Main with File I/O
*/
int main() {
    printHeader("NEWTON'S BACKWARD INTERPOLATION", cout);
    
    string inputFile, outputFile;
    cout << "\nEnter input filename: ";
    cin >> inputFile;
    cout << "Enter output filename: ";
    cin >> outputFile;
    
    // Open input file
    ifstream fin(inputFile);
    if (!fin) {
        cerr << "Error: Cannot open input file '" << inputFile << "'\n";
        return 1;
    }
    
    // Read data points
    int n;
    fin >> n;
    if (n <= 0) {
        cerr << "Error: Invalid number of data points\n";
        return 1;
    }
    
    vector<double> xs(n), ys(n);
    for (int i=0; i<n; i++){
        fin >> xs[i] >> ys[i];
    }
    
    // Read interpolation points
    int m;
    fin >> m;
    if (m <= 0) {
        cerr << "Error: Invalid number of interpolation points\n";
        return 1;
    }
    
    vector<double> xInterpolate(m);
    for (int i=0; i<m; i++){
        fin >> xInterpolate[i];
    }
    
    // Check for additional data point
    double xAdditional, yAdditional;
    bool hasAdditional = false;
    if (fin >> xAdditional >> yAdditional) {
        hasAdditional = true;
    }
    fin.close();
    
    // Verify equal spacing
    if (n > 1) {
        double h = xs[1] - xs[0];
        for (int i=2; i<n; i++){
            if (fabs((xs[i] - xs[i-1]) - h) > 1e-9) {
                cerr << "Error: Data points are not equally spaced!\n";
                return 1;
            }
        }
        cout << "\nNumber of data points: " << n << "\n";
        cout << "Step size (h): " << fixed << setprecision(6) << h << "\n";
    }
    
    // Open output file
    ofstream fout(outputFile);
    if (!fout) {
        cerr << "Error: Cannot create output file '" << outputFile << "'\n";
        return 1;
    }
    
    printHeader("NEWTON'S BACKWARD INTERPOLATION", fout);
    fout << "\nNumber of data points: " << n << "\n";
    if (n > 1) {
        fout << "Step size (h): " << fixed << setprecision(6) << (xs[1] - xs[0]) << "\n";
    }
    
    // Print data points table
    printDataTable(xs, ys, cout);
    printDataTable(xs, ys, fout);
    
    // Build and print backward difference table
    vector<vector<double>> diffTable = buildBackwardDiffTable(xs, ys);
    printBackwardDiffTable(xs, diffTable, cout);
    printBackwardDiffTable(xs, diffTable, fout);
    
    // Perform interpolation
    printHeader("  INTERPOLATION RESULTS", cout);
    printHeader("  INTERPOLATION RESULTS", fout);
    
    vector<double> interpolatedValues(m);
    processInterpolation(xs, diffTable, xInterpolate, interpolatedValues, cout, fout);
    
    cout << "====================================\n";
    fout << "====================================\n";
    
    // Process additional point if provided
    if (hasAdditional) {
        vector<double> xsNew = xs;
        vector<double> ysNew = ys;
        
        // Insert in sorted order
        int insertPos = n;
        for (int i=0; i<n; i++){
            if (xAdditional < xs[i]) {
                insertPos = i;
                break;
            }
        }
        xsNew.insert(xsNew.begin() + insertPos, xAdditional);
        ysNew.insert(ysNew. begin() + insertPos, yAdditional);
        
        // Verify equal spacing
        int nNew = n + 1;
        double hNew = xsNew[1] - xsNew[0];
        bool validSpacing = true;
        for (int i=2; i<nNew; i++){
            if (fabs((xsNew[i] - xsNew[i-1]) - hNew) > 1e-9) {
                validSpacing = false;
                break;
            }
        }
        
        if (validSpacing) {
            printHeader("  WITH ADDITIONAL DATA POINT", cout);
            printHeader("  WITH ADDITIONAL DATA POINT", fout);
            
            cout << "Additional point:  x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            cout << "New number of data points: " << nNew << "\n";
            cout << "New step size (h): " << setprecision(6) << hNew << "\n";
            
            fout << "Additional point: x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            fout << "New number of data points: " << nNew << "\n";
            fout << "New step size (h): " << setprecision(6) << hNew << "\n";
            
            // Print updated data points table
            printDataTable(xsNew, ysNew, cout);
            printDataTable(xsNew, ysNew, fout);
            
            // Build new difference table
            vector<vector<double>> diffTableNew = buildBackwardDiffTable(xsNew, ysNew);
            printBackwardDiffTable(xsNew, diffTableNew, cout);
            printBackwardDiffTable(xsNew, diffTableNew, fout);
            
            // Re-interpolate and show differences
            printHeader("  UPDATED INTERPOLATION RESULTS", cout);
            printHeader("  UPDATED INTERPOLATION RESULTS", fout);
            
            for (int i=0; i<m; i++){
                double resultOld = interpolatedValues[i];
                double resultNew = newtonBackwardWithTable(xsNew, diffTableNew, xInterpolate[i]);
                double absError = fabs(resultNew - resultOld);
                double relError = (fabs(resultNew) > 1e-15) ? (absError / fabs(resultNew)) * 100.0 : 0.0;
                
                cout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                cout << "         Old y = " << setprecision(6) << resultOld << "\n";
                cout << "         New y = " << setprecision(6) << resultNew << "\n";
                cout << "         Absolute Difference: " << scientific << setprecision(6) << absError << "\n";
                cout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
                
                fout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                fout << "         Old y = " << setprecision(6) << resultOld << "\n";
                fout << "         New y = " << setprecision(6) << resultNew << "\n";
                fout << "         Absolute Difference: " << scientific << setprecision(6) << absError << "\n";
                fout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
            }
            
            cout << "====================================\n";
            fout << "====================================\n";
        } else {
            printHeader("WARNING: Additional point breaks equal spacing!", cout);
            cout << "Cannot use Newton's Backward Interpolation.\n";
            cout << "====================================\n";
            
            printHeader("WARNING: Additional point breaks equal spacing!", fout);
            fout << "Cannot use Newton's Backward Interpolation.\n";
            fout << "====================================\n";
        }
    }
    
    fout.close();
    cout << "\nResults have been written to '" << outputFile << "'\n";
    cout << "Program executed successfully!\n";
    
    return 0;
}
```

---

## Backward Interpolation Input
**Input1 (input1.txt):**   
```
5
0.0 1.0
0.5 1.6487
1.0 2.7183
1.5 4.4817
2.0 7.3891
3
1.25
1.75
1.85
```
**Input2 (input2.txt):**   
```
5
0.0 1.0
0.5 1.6487
1.0 2.7183
1.5 4.4817
2.0 7.3891
3
1.25
1.75
1.85
2.5 12.1825
```
---

## Backward Interpolation Output
**Output1 (output1.txt):** 
```

====================================
NEWTON'S BACKWARD INTERPOLATION
====================================

Number of data points: 5
Step size (h): 0.500000

====================================
  DATA POINTS TABLE
====================================
    x	|    0.0000	|    0.5000	|    1.0000	|    1.5000	|    2.0000	|
----------------------------------------------------------------------
    y	|  1.000000	|  1.648700	|  2.718300	|  4.481700	|  7.389100	|
====================================

====================================
  BACKWARD DIFFERENCE TABLE
====================================
           x              y       Nabla^1y       Nabla^2y       Nabla^3y       Nabla^4y
---------------------------------------------------------------------------------------
      0.0000       1.000000
      0.5000       1.648700       0.648700
      1.0000       2.718300       1.069600       0.420900
      1.5000       4.481700       1.763400       0.693800       0.272900
      2.0000       7.389100       2.907400       1.144000       0.450200       0.177300
====================================

====================================
  INTERPOLATION RESULTS
====================================
Point 1: x = 1.250000
         y = 3.489293
Point 2: x = 1.750000
         y = 5.757337
Point 3: x = 1.850000
         y = 6.362852
====================================

```

**Output2 (output2.txt):** 
```

====================================
NEWTON'S BACKWARD INTERPOLATION
====================================

Number of data points: 5
Step size (h): 0.500000

====================================
  DATA POINTS TABLE
====================================
    x	|    0.0000	|    0.5000	|    1.0000	|    1.5000	|    2.0000	|
----------------------------------------------------------------------
    y	|  1.000000	|  1.648700	|  2.718300	|  4.481700	|  7.389100	|
====================================

====================================
  BACKWARD DIFFERENCE TABLE
====================================
           x              y       Nabla^1y       Nabla^2y       Nabla^3y       Nabla^4y
---------------------------------------------------------------------------------------
      0.0000       1.000000
      0.5000       1.648700       0.648700
      1.0000       2.718300       1.069600       0.420900
      1.5000       4.481700       1.763400       0.693800       0.272900
      2.0000       7.389100       2.907400       1.144000       0.450200       0.177300
====================================

====================================
  INTERPOLATION RESULTS
====================================
Point 1: x = 1.250000
         y = 3.489293
Point 2: x = 1.750000
         y = 5.757337
Point 3: x = 1.850000
         y = 6.362852
====================================

====================================
  WITH ADDITIONAL DATA POINT
====================================
Additional point: x = 2.500000, y = 12.182500
New number of data points: 6
New step size (h): 0.500000

====================================
  DATA POINTS TABLE
====================================
    x	|    0.0000	|    0.5000	|    1.0000	|    1.5000	|    2.0000	|    2.5000	|
-----------------------------------------------------------------------------------
    y	|  1.000000	|  1.648700	|  2.718300	|  4.481700	|  7.389100	| 12.182500	|
====================================

====================================
  BACKWARD DIFFERENCE TABLE
====================================
           x              y       Nabla^1y       Nabla^2y       Nabla^3y       Nabla^4y       Nabla^5y
------------------------------------------------------------------------------------------------------
      0.0000       1.000000
      0.5000       1.648700       0.648700
      1.0000       2.718300       1.069600       0.420900
      1.5000       4.481700       1.763400       0.693800       0.272900
      2.0000       7.389100       2.907400       1.144000       0.450200       0.177300
      2.5000      12.182500       4.793400       1.886000       0.742000       0.291800       0.114500
====================================

====================================
  UPDATED INTERPOLATION RESULTS
====================================
Point 1: x = 1.250000
         Old y = 3.489293
         New y = 3.490635
         Absolute Difference: 1.341797e-03
         Relative Difference: 0.0384%

Point 2: x = 1.750000
         Old y = 5.757337
         New y = 5.754206
         Absolute Difference: 3.130859e-03
         Relative Difference: 0.0544%

Point 3: x = 1.850000
         Old y = 6.362852
         New y = 6.359449
         Absolute Difference: 3.402969e-03
         Relative Difference: 0.0535%

====================================
```
---

# Newton's Divided Difference Interpolation
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Interpolation%20and%20Approximation/Newton's%20Divided%20Difference%20Interpolation/)

## Divided Difference Interpolation Theory
### Mathematical Foundation

Given a set of data points $(x_0, y_0), (x_1, y_1), ..., (x_n, y_n)$ with **arbitrary spacing**, Newton's Divided Difference Interpolation formula is:

$$P(x) = f[x_0] + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1) + ... + f[x_0,... ,x_n](x-x_0).. .(x-x_{n-1})$$

where $f[x_0, x_1, ..., x_k]$ denotes the $k$-th order divided difference. 

**Divided Difference Operator:**

The divided differences are computed recursively: 

**Zeroth order (function values):**
$$f[x_i] = y_i$$

**First order:**
$$f[x_i, x_{i+1}] = \frac{f[x_{i+1}] - f[x_i]}{x_{i+1} - x_i}$$

**Second order:**
$$f[x_i, x_{i+1}, x_{i+2}] = \frac{f[x_{i+1}, x_{i+2}] - f[x_i, x_{i+1}]}{x_{i+2} - x_i}$$

**General k-th order:**
$$f[x_i, x_{i+1}, ..., x_{i+k}] = \frac{f[x_{i+1}, ..., x_{i+k}] - f[x_i, ..., x_{i+k-1}]}{x_{i+k} - x_i}$$

**Divided Difference Table:**

| $x$ | $f(x)$ | $f[x_0,x_1]$ | $f[x_0,x_1,x_2]$ | $f[x_0,x_1,x_2,x_3]$ | $f[x_0,... ,x_4]$ |
|-----|--------|-------------|-----------------|---------------------|------------------|
| $x_0$ | $f[x_0]$ | $f[x_0,x_1]$ | $f[x_0,x_1,x_2]$ | $f[x_0,x_1,x_2,x_3]$ | $f[x_0,... ,x_4]$ |
| $x_1$ | $f[x_1]$ | $f[x_1,x_2]$ | $f[x_1,x_2,x_3]$ | $f[x_1,x_2,x_3,x_4]$ | |
| $x_2$ | $f[x_2]$ | $f[x_2,x_3]$ | $f[x_2,x_3,x_4]$ | | |
| $x_3$ | $f[x_3]$ | $f[x_3,x_4]$ | | | |
| $x_4$ | $f[x_4]$ | | | | |

### Algorithm Steps

1. **Read Data Points**: No spacing validation needed (works with any spacing)
2. **Build Divided Difference Table**:
   - Initialize first column with $f(x)$ values
   - Compute successive divided differences using:  $f[x_i,... ,x_{i+k}] = \frac{f[x_{i+1},...,x_{i+k}] - f[x_i,... ,x_{i+k-1}]}{x_{i+k} - x_i}$
   - Check for duplicate x-values (would cause division by zero)
3. **Apply Newton's Divided Difference Formula**:
   - Start with $P(x) = f[x_0]$
   - Add terms: $f[x_0,... ,x_j] \cdot (x-x_0) \cdot ...  \cdot (x-x_{j-1})$ for $j = 1, 2, ..., n$
4. **Return Interpolated Value**

### Complexity Analysis

- **Time Complexity**:
  - Building divided difference table: $O(n^2)$ where $n$ is the number of data points
  - Single interpolation: $O(n)$
  - Total for $m$ interpolations: $O(n^2 + mn)$
  
- **Space Complexity**:  $O(n^2)$ for storing the divided difference table

- **Accuracy**: 
  - Provides exact results for polynomial data of degree $\leq n-1$
  - Works with **any spacing** between data points
  - Uniform accuracy across the entire data range
  - More versatile than Forward/Backward interpolation
---

## Divided Difference Interpolation Code
```cpp
#include <bits/stdc++.h>
using namespace std;

/*
   Print Data Points Table
*/
void printDataTable(const vector<double>& xs, const vector<double>& ys, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  DATA POINTS TABLE\n";
    out << "====================================\n";
    
    // Print x values
    out << setw(10) << "x";
    for (int i=0; i<n; i++){
        out << setw(15) << fixed << setprecision(4) << xs[i];
    }
    out << "\n";
    
    // Print separator
    out << string(10 + 15*n, '-') << "\n";
    
    // Print y values
    out << setw(10) << "y";
    for (int i=0; i<n; i++){
        out << setw(15) << fixed << setprecision(6) << ys[i];
    }
    out << "\n";
    out << "====================================\n";
}

/*
   Build Divided Difference Table
*/
vector<vector<double>> buildDividedDiffTable(const vector<double>& xs, const vector<double>& ys) {
    int n = (int)xs.size();
    vector<vector<double>> diff(n, vector<double>(n));
    
    // First column is y values
    for (int i=0; i<n; i++) diff[i][0] = ys[i];
    
    // Build divided difference table
    for (int j=1; j<n; j++){
        for (int i=0; i<n-j; i++){
            double denom = xs[i+j] - xs[i];
            if (fabs(denom) < 1e-15) {
                throw runtime_error("Duplicate x-values encountered");
            }
            diff[i][j] = (diff[i+1][j-1] - diff[i][j-1]) / denom;
        }
    }
    return diff;
}

/*
   Print Divided Difference Table
*/
void printDividedDiffTable(const vector<double>& xs, const vector<vector<double>>& diff, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  DIVIDED DIFFERENCE TABLE\n";
    out << "====================================\n";
    out << setw(12) << "x" << setw(15) << "f(x)";
    for (int j=1; j<n; j++){
        out << setw(15) << ("f[x" + to_string(0) + "..." + "x" + to_string(j) + "]");
    }
    out << "\n" << string(12 + 15*n, '-') << "\n";
    
    for (int i=0; i<n; i++){
        out << fixed << setprecision(4) << setw(12) << xs[i];
        for (int j=0; j<n-i; j++){
            out << setprecision(6) << setw(15) << diff[i][j];
        }
        out << "\n";
    }
    out << "====================================\n";
}

/*
   Newton's Divided Difference Interpolation using pre-built table
*/
double newtonDividedDifferenceWithTable(const vector<double>& xs, const vector<vector<double>>& diff, double x) {
    int n = (int)xs.size();
    if (n == 1) return diff[0][0];
    
    double result = diff[0][0];
    double term = 1.0;
    
    for (int j=1; j<n; j++){
        term *= (x - xs[j-1]);
        result += diff[0][j] * term;
    }
    return result;
}

/*
   Process and output interpolation results
*/
void processInterpolation(const vector<double>& xs, const vector<vector<double>>& diff,
                         const vector<double>& xInterpolate, vector<double>& results,
                         ostream& cout_stream, ostream& fout) {
    int m = xInterpolate.size();
    int n = xs.size();
    
    for (int i=0; i<m; i++){
        double result = newtonDividedDifferenceWithTable(xs, diff, xInterpolate[i]);
        results[i] = result;
        
        bool isExtrap = (xInterpolate[i] < xs[0] || xInterpolate[i] > xs[n-1]);
        string extrapNote = isExtrap ? " (Extrapolation)" : "";
        
        string output = "Point " + to_string(i+1) + ": x = ";
        cout_stream << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        cout_stream << "         y = " << setprecision(6) << result << "\n";
        
        fout << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        fout << "         y = " << setprecision(6) << result << "\n";
    }
}

/*
   Print section header
*/
void printHeader(const string& title, ostream& out) {
    out << "\n====================================\n";
    out << title << "\n";
    out << "====================================\n";
}

/*
   Main with File I/O
*/
int main() {
    printHeader("NEWTON'S DIVIDED DIFFERENCE INTERPOLATION", cout);
    
    string inputFile, outputFile;
    cout << "\nEnter input filename:  ";
    cin >> inputFile;
    cout << "Enter output filename: ";
    cin >> outputFile;
    
    // Open input file
    ifstream fin(inputFile);
    if (!fin) {
        cerr << "Error: Cannot open input file '" << inputFile << "'\n";
        return 1;
    }
    
    // Read data points
    int n;
    fin >> n;
    if (n <= 0) {
        cerr << "Error: Invalid number of data points\n";
        return 1;
    }
    
    vector<double> xs(n), ys(n);
    for (int i=0; i<n; i++){
        fin >> xs[i] >> ys[i];
    }
    
    // Read interpolation points
    int m;
    fin >> m;
    if (m <= 0) {
        cerr << "Error: Invalid number of interpolation points\n";
        return 1;
    }
    
    vector<double> xInterpolate(m);
    for (int i=0; i<m; i++){
        fin >> xInterpolate[i];
    }
    
    // Check for additional data point
    double xAdditional, yAdditional;
    bool hasAdditional = false;
    if (fin >> xAdditional >> yAdditional) {
        hasAdditional = true;
    }
    fin.close();
    
    cout << "\nNumber of data points: " << n << "\n";
    
    // Open output file
    ofstream fout(outputFile);
    if (!fout) {
        cerr << "Error: Cannot create output file '" << outputFile << "'\n";
        return 1;
    }
    
    printHeader("NEWTON'S DIVIDED DIFFERENCE INTERPOLATION", fout);
    fout << "\nNumber of data points: " << n << "\n";
    
    // Print data points table
    printDataTable(xs, ys, cout);
    printDataTable(xs, ys, fout);
    
    // Build and print divided difference table
    try {
        vector<vector<double>> diffTable = buildDividedDiffTable(xs, ys);
        printDividedDiffTable(xs, diffTable, cout);
        printDividedDiffTable(xs, diffTable, fout);
        
        // Perform interpolation
        printHeader("  INTERPOLATION RESULTS", cout);
        printHeader("  INTERPOLATION RESULTS", fout);
        
        vector<double> interpolatedValues(m);
        processInterpolation(xs, diffTable, xInterpolate, interpolatedValues, cout, fout);
        
        cout << "====================================\n";
        fout << "====================================\n";
        
        // Process additional point if provided
        if (hasAdditional) {
            vector<double> xsNew = xs;
            vector<double> ysNew = ys;
            
            // Insert in sorted order
            int insertPos = n;
            for (int i=0; i<n; i++){
                if (xAdditional < xs[i]) {
                    insertPos = i;
                    break;
                }
            }
            xsNew.insert(xsNew.begin() + insertPos, xAdditional);
            ysNew.insert(ysNew. begin() + insertPos, yAdditional);
            
            int nNew = n + 1;
            
            printHeader("  WITH ADDITIONAL DATA POINT", cout);
            printHeader("  WITH ADDITIONAL DATA POINT", fout);
            
            cout << "Additional point: x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            cout << "New number of data points: " << nNew << "\n";
            
            fout << "Additional point: x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            fout << "New number of data points: " << nNew << "\n";
            
            // Print updated data points table
            printDataTable(xsNew, ysNew, cout);
            printDataTable(xsNew, ysNew, fout);
            
            // Build new divided difference table
            vector<vector<double>> diffTableNew = buildDividedDiffTable(xsNew, ysNew);
            printDividedDiffTable(xsNew, diffTableNew, cout);
            printDividedDiffTable(xsNew, diffTableNew, fout);
            
            // Re-interpolate and show differences
            printHeader("  UPDATED INTERPOLATION RESULTS", cout);
            printHeader("  UPDATED INTERPOLATION RESULTS", fout);
            
            for (int i=0; i<m; i++){
                double resultOld = interpolatedValues[i];
                double resultNew = newtonDividedDifferenceWithTable(xsNew, diffTableNew, xInterpolate[i]);
                double absError = fabs(resultNew - resultOld);
                double relError = (fabs(resultNew) > 1e-15) ? (absError / fabs(resultNew)) * 100.0 : 0.0;
                
                cout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                cout << "         Old y = " << setprecision(6) << resultOld << "\n";
                cout << "         New y = " << setprecision(6) << resultNew << "\n";
                cout << "         Absolute Difference: " << scientific << setprecision(6) << absError << "\n";
                cout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
                
                fout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                fout << "         Old y = " << setprecision(6) << resultOld << "\n";
                fout << "         New y = " << setprecision(6) << resultNew << "\n";
                fout << "         Absolute Difference: " << scientific << setprecision(6) << absError << "\n";
                fout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
            }
            
            cout << "====================================\n";
            fout << "====================================\n";
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        fout << "Error: " << e.what() << "\n";
        return 1;
    }
    
    fout.close();
    cout << "\nResults have been written to '" << outputFile << "'\n";
    cout << "Program executed successfully!\n";
    
    return 0;
}
```

---

## Divided Difference Interpolation Input
**Input1 (input1.txt):**   
```
5
0.0 1.0
0.7 2.014
1.3 3.669
2.0 7.389
2.8 16.445
3
0.5
1.5
2.5
```
**Input2 (input2.txt):**   
```
5
0.0 1.0
0.7 2.014
1.3 3.669
2.0 7.389
2.8 16.445
3
0.5
1.5
2.5
3.5 33.115
```
---

## Divided Difference Interpolation Output
**Output1 (output1.txt):** 
```

====================================
NEWTON'S DIVIDED DIFFERENCE INTERPOLATION
====================================

Number of data points: 5

====================================
  DATA POINTS TABLE
====================================
         x         0.0000         0.7000         1.3000         2.0000         2.8000
-------------------------------------------------------------------------------------
         y       1.000000       2.014000       3.669000       7.389000      16.445000
====================================

====================================
  DIVIDED DIFFERENCE TABLE
====================================
           x           f(x)     f[x0...x1]     f[x0...x2]     f[x0...x3]     f[x0...x4]
---------------------------------------------------------------------------------------
      0.0000       1.000000       1.448571       1.007509       0.479304       0.175366
      0.7000       2.014000       2.758333       1.966117       0.970330
      1.3000       3.669000       5.314286       4.003810
      2.0000       7.389000      11.320000
      2.8000      16.445000
====================================

====================================
  INTERPOLATION RESULTS
====================================
Point 1: x = 0.500000
         y = 1.640835
Point 2: x = 1.500000
         y = 4.475857
Point 3: x = 2.500000
         y = 12.216951
====================================

```

**Output2 (output2.txt):** 
```

====================================
NEWTON'S DIVIDED DIFFERENCE INTERPOLATION
====================================

Number of data points: 5

====================================
  DATA POINTS TABLE
====================================
         x         0.0000         0.7000         1.3000         2.0000         2.8000
-------------------------------------------------------------------------------------
         y       1.000000       2.014000       3.669000       7.389000      16.445000
====================================

====================================
  DIVIDED DIFFERENCE TABLE
====================================
           x           f(x)     f[x0...x1]     f[x0...x2]     f[x0...x3]     f[x0...x4]
---------------------------------------------------------------------------------------
      0.0000       1.000000       1.448571       1.007509       0.479304       0.175366
      0.7000       2.014000       2.758333       1.966117       0.970330
      1.3000       3.669000       5.314286       4.003810
      2.0000       7.389000      11.320000
      2.8000      16.445000
====================================

====================================
  INTERPOLATION RESULTS
====================================
Point 1: x = 0.500000
         y = 1.640835
Point 2: x = 1.500000
         y = 4.475857
Point 3: x = 2.500000
         y = 12.216951
====================================

====================================
  WITH ADDITIONAL DATA POINT
====================================
Additional point: x = 3.500000, y = 33.115000
New number of data points: 6

====================================
  DATA POINTS TABLE
====================================
         x         0.0000         0.7000         1.3000         2.0000         2.8000         3.5000
----------------------------------------------------------------------------------------------------
         y       1.000000       2.014000       3.669000       7.389000      16.445000      33.115000
====================================

====================================
  DIVIDED DIFFERENCE TABLE
====================================
           x           f(x)     f[x0...x1]     f[x0...x2]     f[x0...x3]     f[x0...x4]     f[x0...x5]
------------------------------------------------------------------------------------------------------
      0.0000       1.000000       1.448571       1.007509       0.479304       0.175366       0.051518
      0.7000       2.014000       2.758333       1.966117       0.970330       0.355680
      1.3000       3.669000       5.314286       4.003810       1.966234
      2.0000       7.389000      11.320000       8.329524
      2.8000      16.445000      23.814286
      3.5000      33.115000
====================================

====================================
  UPDATED INTERPOLATION RESULTS
====================================
Point 1: x = 0.500000
         Old y = 1.640835
         New y = 1.655054
         Absolute Difference: 1.421903e-02
         Relative Difference: 0.8591%

Point 2: x = 1.500000
         Old y = 4.475857
         New y = 4.483894
         Absolute Difference: 8.036841e-03
         Relative Difference: 0.1792%

Point 3: x = 2.500000
         Old y = 12.216951
         New y = 12.175221
         Absolute Difference: 4.172975e-02
         Relative Difference: 0.3427%

====================================
```
---
---

# Numerical Integration
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-red?style=for-the-badge)](./Numerical%20Integration/)
## 📐 Mathematical Foundation
Numerical Integration (Numerical Quadrature) is used to approximate definite integrals when:

an analytical integral is difficult or impossible,
the function values are given in tabular form,
or a quick approximation is required.
Simpson’s rules are popular Newton–Cotes methods that assume equally spaced x-values and approximate the integrand using low-degree polynomials.
---

# Simpson's One-third Rule
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Numerical%20Integration/Simpson's%20One-third%20Rule/)

## Simpson One-third Theory
### 📌 Mathematical Formula

Divide the interval \([a, b]\) into **n** equal subintervals (**n** must be even).

- **Step size:**  
  `h = (b - a) / n`
- **Points:**  
  `x_i = a + i * h`, where `i = 0, 1, ..., n`

Then:

```
∫[a to b] f(x) dx ≈ (h/3) × [ y₀ + yₙ + 4(y₁ + y₃ + ... + yₙ₋₁) + 2(y₂ + y₄ + ... + yₙ₋₂) ]
```
where `y_i = f(x_i)`.


### ✅ Validity Condition

- `n` must be **even**
- data must be **equally spaced**


### 🧾 Algorithm Steps

1. Read `n`
2. Read `a` and `b`
3. Read `y0..yn` (total `n+1` values)
4. Check `n` is even (otherwise print error for that case)
5. Compute `h = (b-a)/n`
6. Compute:
   - `sum_odd = y1 + y3 + ... + y(n-1)`
   - `sum_even = y2 + y4 + ... + y(n-2)`
7. Compute integral using Simpson’s 1/3 formula
8. Print result in a clear formatted way
---

## Simpson One-third Code
```cpp
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

int main() {
    int n;
    double a, b;
    cin >> n >> a >> b;

    vector<double> y(n+1);
    for(int i=0;i<=n;i++) cin >> y[i];

    if(n % 2 != 0){
        cout << "n must be even";
        return 0;
    }

    double h = (b-a)/n;
    double odd=0, even=0;

    for(int i=1;i<n;i++){
        if(i%2==0) even+=y[i];
        else odd+=y[i];
    }

    double result = (h/3)*(y[0]+y[n]+4*odd+2*even);
    cout << fixed << setprecision(6) << result;
    return 0;
}

```

---

## Simpson One-third Input
**Input (input.txt):**   
```
6
0 1
1 0.6944 0.4444 0.25 0.1111 0
```

## Simpson One-third Output
**Output (output.txt):** 
```
0.333333
```

---

# Simpson's Three-eighths Rule
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Numerical%20Integration/Simpson's%20Three-eighths%20Rule/)

## Simpson Three-eighths Theory
### 📌 Mathematical Formula

> Divide the interval ([a, b]) into **n** equal subintervals (**n** must be a multiple of 3).

- **Step size:**  
  `h = (b - a) / n`
- **Points:**  
  `x_i = a + i * h`, where `i = 0, 1, ..., n`

Then:

```
∫[a to b] f(x) dx ≈ (3h / 8) * [
    y₀ + yₙ
    + 3(y₁ + y₂ + y₄ + y₅ + ...)
    + 2(y₃ + y₆ + y₉ + ...)
]
```

where `yᵢ = f(xᵢ)`.


### ✅ Validity Condition

- `n` must be a **multiple of 3**
- data must be **equally spaced**


### 🧾 Algorithm Steps

1. Read `n`
2. Read `a` and `b`
3. Read `y0..yn` (total `n+1` values)
4. Check `n` is divisible by 3 (otherwise print error for that case)
5. Compute `h = (b - a) / n`
6. Compute:
   - `sum_3 = y1 + y2 + y4 + y5 + ...` (indices not divisible by 3)
   - `sum_2 = y3 + y6 + y9 + ...` (indices divisible by 3, excluding endpoints)
7. Compute integral using Simpson’s 3/8 formula
8. Print result in a clear formatted way
---

## Simpson Three-eighths Code
```cpp
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

int main() {
    int n;
    double a, b;
    cin >> n >> a >> b;

    vector<double> y(n+1);
    for(int i=0;i<=n;i++) cin >> y[i];

    if(n % 3 != 0){
        cout << "n must be multiple of 3";
        return 0;
    }

    double h = (b-a)/n;
    double sum3=0, sum2=0;

    for(int i=1;i<n;i++){
        if(i%3==0) sum2+=y[i];
        else sum3+=y[i];
    }

    double result = (3*h/8)*(y[0]+y[n]+3*sum3+2*sum2);
    cout << fixed << setprecision(6) << result;
    return 0;
}

```

---

## Simpson Three-eighths Input
**Input (input.txt):**   
```
6
0 1
1 0.6944 0.4444 0.25 0.1111 0
```
---

## Simpson Three-eighths Output
**Output (output.txt):** 
```
0.333333
```
---
---

# Numerical Differentiation
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-red?style=for-the-badge)](./Numerical%20Differentiation/)
## 📐 Mathematical Foundation
### What is Numerical Differentiation?

**Numerical differentiation** is the process of estimating the derivative of a function using values at discrete points. When $f(x)$ is known only at certain $x$ values, we use difference quotients to approximate $f'(x)$ and $f''(x)$.

**Formal Definition**: The derivative at $x$ is defined as:

$$
f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}
$$

In practice, with finite $h$, we use difference formulas.

### Finite Difference Approximations

#### **Forward Difference**

$$
f'(x_0) \approx \frac{f(x_1) - f(x_0)}{h}
$$

#### **Backward Difference**

$$
f'(x_n) \approx \frac{f(x_n) - f(x_{n-1})}{h}
$$

#### **Central Difference**

$$
f'(x_i) \approx \frac{f(x_{i+1}) - f(x_{i-1})}{2h}
$$

**Key Concepts**:
- $h$ is the spacing between $x$ values (assumed equal for standard formulas)
- Higher-order differences can be used for better accuracy
- Formulas can be extended for non-uniform spacing (see advanced texts)

### Error Analysis

Numerical differentiation introduces **truncation error** (from Taylor expansion) and **round-off error** (from finite precision arithmetic).

**Error for Forward/Backward Difference**:
$$
f'(x) = \frac{f(x+h) - f(x)}{h} - \frac{h}{2}f''(\xi)
$$
for some $\xi$ in $[x, x+h]$. Error is $O(h)$.

**Error for Central Difference**:
$$
f'(x) = \frac{f(x+h) - f(x-h)}{2h} - \frac{h^2}{6}f'''(\xi)
$$
for some $\xi$ in $[x-h, x+h]$. Error is $O(h^2)$.

**Key Insights**:
- Smaller $h$ reduces truncation error but increases round-off error
- Central difference is generally more accurate
- For noisy data, differentiation amplifies noise—smoothing may be needed
---

# First and Second Order Derivative based on Forward Interpolation
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Numerical%20Differentiation/First%20and%20Second%20Order%20Derivative%20based%20on%20Forward%20Interpolation/)

## Forward Interpolation Derivative Theory
### Mathematical Foundation

Given $n$ equally spaced data points $(x_0, y_0), (x_1, y_1), ..., (x_{n-1}, y_{n-1})$, the forward difference table is constructed. The first and second derivatives at a point $x_p$ (close to $x_0$) are approximated as:

- **First Derivative:**
   $$
   f'(x_p) \approx \frac{1}{h} \left[\Delta y_0 + \frac{2u-1}{2} \Delta^2 y_0 + \frac{3u^2-6u+2}{6} \Delta^3 y_0 + ...\right]
   $$
- **Second Derivative:**
   $$
   f''(x_p) \approx \frac{1}{h^2} \left[\Delta^2 y_0 + (u-1)\Delta^3 y_0 + \frac{6u^2-18u+11}{12} \Delta^4 y_0 + ...\right]
   $$
where $u = \frac{x_p - x_0}{h}$ and $h$ is the spacing between $x$ values.

#### Forward Difference Table

| $x$ | $y$ | $\Delta y$ | $\Delta^2 y$ | $\Delta^3 y$ | $\Delta^4 y$ |
|-----|-----|-----------|-------------|-------------|-------------|
| $x_0$ | $y_0$ | $\Delta y_0$ | $\Delta^2 y_0$ | $\Delta^3 y_0$ | $\Delta^4 y_0$ |
| $x_1$ | $y_1$ | $\Delta y_1$ | $\Delta^2 y_1$ | $\Delta^3 y_1$ | |
| $x_2$ | $y_2$ | $\Delta y_2$ | $\Delta^2 y_2$ | | |
| $x_3$ | $y_3$ | $\Delta y_3$ | | | |
| $x_4$ | $y_4$ | | | | |

### Algorithm Steps

1. **Read** $n$, $x_i$, $y_i$, and $x_p$ from input file.
2. **Build** the forward difference table for $y$ values.
3. **Compute** $u = (x_p - x_0)/h$.
4. **Calculate** the first and second derivatives using the above formulas.
5. **Output** results in a formatted table to `output.txt`.

### Complexity Analysis

- **Time Complexity:**
   - Building difference table: $O(n^2)$
   - Derivative calculation: $O(n)$
- **Space Complexity:** $O(n^2)$

---

## Forward Interpolation Derivative Code
```cpp
#include <bits/stdc++.h>
using namespace std;

// Utility: factorial
static double factorial(int n)
{
    double f = 1.0;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

// Utility: Forward Difference Table
vector<vector<double>> forwardDiff(const vector<double> &y)
{
    int n = y.size();
    vector<vector<double>> d(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        d[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            d[i][j] = d[i + 1][j - 1] - d[i][j - 1];

    return d;
}

// Forward First Derivative
static double forwardFirstDerivative(const vector<double> &x,
                                     const vector<double> &y,
                                     double xp)
{
    double h = x[1] - x[0];
    auto d = forwardDiff(y);
    double u = (xp - x[0]) / h;
    double res = d[0][1];

    if (x.size() >= 3)
        res += ((2 * u - 1) / 2.0) * d[0][2];
    if (x.size() >= 4)
        res += ((3 * u * u - 6 * u + 2) / 6.0) * d[0][3];

    return res / h;
}

// Forward Second Derivative
static double forwardSecondDerivative(const vector<double> &x,
                                      const vector<double> &y,
                                      double xp)
{
    double h = x[1] - x[0];
    auto d = forwardDiff(y);
    double u = (xp - x[0]) / h;
    double res = d[0][2];

    if (x.size() >= 4)
        res += (u - 1) * d[0][3];
    if (x.size() >= 5)
        res += ((6 * u * u - 18 * u + 11) / 12.0) * d[0][4];
        
    return res / (h * h);
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin.is_open())
    {
        cerr << "Error: Could not open input.txt" << endl;
        return 1;
    }
    if (!fout.is_open())
    {
        cerr << "Error: Could not open output.txt" << endl;
        return 1;
    }

    int n;
    fin >> n;

    vector<double> x(n), y(n);
    for (int i = 0; i < n; ++i)
        fin >> x[i];
    for (int i = 0; i < n; ++i)
        fin >> y[i];

    double xp;
    fin >> xp;

    fout << fixed << setprecision(6);
    fout << "Numerical Differentiation using Forward Interpolation\n";
    fout << "---------------------------------------------------\n";
    
    fout << "Given data points (x, y):\n";
    for (int i = 0; i < n; ++i)
        fout << "  (" << setw(8) << x[i] << ", " << setw(10) << y[i] << ")\n";
    fout << "\n";
    
    fout << "Point of differentiation (xp): " << xp << "\n\n";

    double first = forwardFirstDerivative(x, y, xp);
    double second = forwardSecondDerivative(x, y, xp);

    fout << "First Derivative at x = " << xp << " : " << first << "\n";
    fout << "Second Derivative at x = " << xp << " : " << second << "\n";
    fout << "\n";

    fin.close();
    fout.close();

    return 0;
}
```

---

## Forward Interpolation Derivative Input
**Input (input.txt):**   
```
5
1 2 3 4 5
1 8 27 64 125
2.5

```
---

## Forward Interpolation Derivative Output
**Output (output.txt):** 
```
Numerical Differentiation using Forward Interpolation
---------------------------------------------------
Given data points (x, y):
  (1.000000,   1.000000)
  (2.000000,   8.000000)
  (3.000000,  27.000000)
  (4.000000,  64.000000)
  (5.000000, 125.000000)

Point of differentiation (xp): 2.500000

First Derivative at x = 2.500000 : 18.750000
Second Derivative at x = 2.500000 : 15.000000
```
---

# First and Second Order Derivative based on Backward Interpolation
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Numerical%20Differentiation/First%20and%20Second%20Order%20Derivative%20based%20on%20Backward%20Interpolation/)

## Backward Interpolation Derivative Theory
### Mathematical Foundation

Given $n$ equally spaced data points $(x_0, y_0), (x_1, y_1), ..., (x_{n-1}, y_{n-1})$, the backward difference table is constructed. The first and second derivatives at a point $x_p$ (close to $x_{n-1}$) are approximated as:

- **First Derivative:**
  $$
  f'(x_p) \approx \frac{1}{h} \left[\nabla y_{n-1} + \frac{2v+1}{2} \nabla^2 y_{n-1} + \frac{3v^2+6v+2}{6} \nabla^3 y_{n-1} + ...\right]
  $$
- **Second Derivative:**
  $$
  f''(x_p) \approx \frac{1}{h^2} \left[\nabla^2 y_{n-1} + (v+1)\nabla^3 y_{n-1} + \frac{6v^2+18v+11}{12} \nabla^4 y_{n-1} + ...\right]
  $$
where $v = \frac{x_p - x_{n-1}}{h}$ and $h$ is the spacing between $x$ values.

#### Backward Difference Table

| $x$ | $y$ | $\nabla y$ | $\nabla^2 y$ | $\nabla^3 y$ | $\nabla^4 y$ |
|-----|-----|-----------|-------------|-------------|-------------|
| $x_0$ | $y_0$ | | | | |
| $x_1$ | $y_1$ | $\nabla y_1$ | | | |
| $x_2$ | $y_2$ | $\nabla y_2$ | $\nabla^2 y_2$ | | |
| $x_3$ | $y_3$ | $\nabla y_3$ | $\nabla^2 y_3$ | $\nabla^3 y_3$ | |
| $x_4$ | $y_4$ | $\nabla y_4$ | $\nabla^2 y_4$ | $\nabla^3 y_4$ | $\nabla^4 y_4$ |

### Algorithm Steps

1. **Read** $n$, $x_i$, $y_i$, and $x_p$ from input file.
2. **Build** the backward difference table for $y$ values.
3. **Compute** $v = (x_p - x_{n-1})/h$.
4. **Calculate** the first and second derivatives using the above formulas.
5. **Output** results in a formatted table to `output.txt`.

### Complexity Analysis

- **Time Complexity:**
  - Building difference table: $O(n^2)$
  - Derivative calculation: $O(n)$
- **Space Complexity:** $O(n^2)$

---

## Backward Interpolation Derivative Code
```cpp
#include <bits/stdc++.h>
using namespace std;

// Utility: factorial
static double factorial(int n)
{
    double f = 1.0;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

// Utility: Backward Difference Table
vector<vector<double>> backwardDiff(const vector<double> &y)
{
    int n = y.size();

    vector<vector<double>> b(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        b[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = n - 1; i >= j; i--)
            b[i][j] = b[i][j - 1] - b[i - 1][j - 1];

    return b;
}

// Backward First Derivative
static double backwardFirstDerivative(const vector<double> &x,
                                      const vector<double> &y,
                                      double xp)
{
    int n = x.size();

    double h = x[1] - x[0];
    auto b = backwardDiff(y);
    double v = (xp - x[n - 1]) / h;
    double res = b[n - 1][1];

    if (n >= 3)
        res += ((2 * v + 1) / 2.0) * b[n - 1][2];
    if (n >= 4)
        res += ((3 * v * v + 6 * v + 2) / 6.0) * b[n - 1][3];

    return res / h;
}

// Backward Second Derivative
static double backwardSecondDerivative(const vector<double> &x,
                                       const vector<double> &y,
                                       double xp)
{
    int n = x.size();

    double h = x[1] - x[0];
    auto b = backwardDiff(y);
    double v = (xp - x[n - 1]) / h;
    double res = b[n - 1][2];

    if (n >= 4)
        res += (v + 1) * b[n - 1][3];
    if (n >= 5)
        res += ((6 * v * v + 18 * v + 11) / 12.0) * b[n - 1][4];
        
    return res / (h * h);
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin.is_open())
    {
        cerr << "Error: Could not open input.txt" << endl;
        return 1;
    }
    if (!fout.is_open())
    {
        cerr << "Error: Could not open output.txt" << endl;
        return 1;
    }

    int n;
    fin >> n;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; ++i)
        fin >> x[i];
    for (int i = 0; i < n; ++i)
        fin >> y[i];

    double xp;
    fin >> xp;

    fout << fixed << setprecision(6);
    fout << "Numerical Differentiation using Backward Interpolation\n";
    fout << "---------------------------------------------------\n";
    
    fout << "Given data points (x, y):\n";
    for (int i = 0; i < n; ++i)
        fout << "  (" << setw(8) << x[i] << ", " << setw(10) << y[i] << ")\n";
    fout << "\n";

    fout << "Point of differentiation (xp): " << xp << "\n\n";

    double first = backwardFirstDerivative(x, y, xp);
    double second = backwardSecondDerivative(x, y, xp);
    fout << "First Derivative at x = " << xp << " : " << first << "\n";
    fout << "Second Derivative at x = " << xp << " : " << second << "\n";
    fout << "\n";

    fin.close();
    fout.close();
    
    return 0;
}

```

---

## Backward Interpolation Derivative Input
**Input (input.txt):**   
```
5
1 2 3 4 5
1 8 27 64 125
4.5
```

---

## Backward Interpolation Derivative Output
**Output (output.txt):** 
```
Numerical Differentiation using Backward Interpolation
---------------------------------------------------
Given data points (x, y):
  (1.000000,   1.000000)
  (2.000000,   8.000000)
  (3.000000,  27.000000)
  (4.000000,  64.000000)
  (5.000000, 125.000000)

Point of differentiation (xp): 4.500000

First Derivative at x = 4.500000 : 60.750000
Second Derivative at x = 4.500000 : 27.000000
```
---
---

# Solution of Differential Equations
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-red?style=for-the-badge)](./Solution%20of%20Differential%20Equations/)
## 📐 Mathematical Foundation
### Ordinary Differential Equations (ODEs)

An ODE relates a function to its derivatives:
```
y'(x) = f(x, y),    y(x_0) = y_0
```
- **Initial Value Problem (IVP):** Given $y(x_0) = y_0$, find $y(x)$ for $x > x_0$.

### Discretization and Step Size

- Choose a step size $h$.
- Define grid points: $x_n = x_0 + n h$ for $n=0,1,2,...,N$.


### Iterative Solution

RK4 advances the solution from $y_n$ to $y_{n+1}$ as follows:
```
k₁ = f(x_n, y_n)
k₂ = f(x_n + h/2, y_n + h·k₁/2)
k₃ = f(x_n + h/2, y_n + h·k₂/2)
k₄ = f(x_n + h,   y_n + h·k₃)
y_{n+1} = y_n + (h/6)·(k₁ + 2·k₂ + 2·k₃ + k₄)
```

### Key Mathematical Concepts

- **Initial Value Problem formulation**
- **First-order ODEs**
- **Discretization of the domain**
- **Incremental/iterative computation**
---

# Runge-Kutta Method (RK4)
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./)

## RK4 Theory
### Initial Value Problem (IVP)

Given:
- ODE: `y' = f(x, y)`
- Initial condition: `y(x₀) = y₀`
- Step size: `h`

We compute `y` at `x₀ + h`, `x₀ + 2h`, … up to the target point.

### RK4 Formula

For each step:

```
k₁ = h * f(xₙ,      yₙ)
k₂ = h * f(xₙ+h/2,  yₙ+k₁/2)
k₃ = h * f(xₙ+h/2,  yₙ+k₂/2)
k₄ = h * f(xₙ+h,    yₙ+k₃)
```

Update:

```
yₙ₊₁ = yₙ + (k₁ + 2k₂ + 2k₃ + k₄) / 6
xₙ₊₁ = xₙ + h
```

### 🧩 Algorithm

1. Read `x₀`, `y₀`
2. Read target `xₙ` and step size `h`
3. Repeat while `x < xₙ`:
   - compute `k₁`, `k₂`, `k₃`, `k₄`
   - update `y` and `x`
4. Print a table of steps and the final answer


### ⏱️ Complexity

If total steps = `N`:

- Time: **O(N)**
- Space: **O(1)** (only a few variables needed)

---

## RK4 Code
```cpp
#include <bits/stdc++.h>
using namespace std;


static double f(int option, double x, double y) {
    switch (option) {
        case 1: return x + y;                 // y' = x + y
        case 2: return x - y;                 // y' = x - y
        case 3: return y - x*x + 1.0;         // y' = y - x^2 + 1
        case 4: return x * y;                 // y' = x*y
        default: return x + y;                // fallback
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int option;
    double x0, y0, xt, h;

    cin >> option;
    cin >> x0 >> y0;
    cin >> xt >> h;

    if (!cin || h <= 0) {
        cout << "Invalid input.\n";
        return 0;
    }
    if (xt < x0) {
        cout << "Target x must be >= x0.\n";
        return 0;
    }

    cout << fixed << setprecision(6);

    cout << "Runge-Kutta 4th Order (RK4)\n";
    cout << "Chosen function option = " << option << "\n";
    cout << "Initial condition: x0 = " << x0 << ", y0 = " << y0 << "\n";
    cout << "Target: x = " << xt << ", step size h = " << h << "\n\n";

    cout << left
         << setw(6)  << "Step"
         << setw(12) << "x"
         << setw(14) << "y"
         << setw(14) << "k1"
         << setw(14) << "k2"
         << setw(14) << "k3"
         << setw(14) << "k4"
         << "\n";

    cout << string(88, '-') << "\n";

    double x = x0, y = y0;
    int step = 0;


    int N = (int)round((xt - x0) / h);
    double x_end = x0 + N * h;

    if (fabs(x_end - xt) > 1e-9) {
        cout << "\nWarning: (xt - x0)/h is not an integer. "
             << "Program will stop at x = " << x_end << " (closest grid point).\n\n";
    }

    for (int i = 0; i < N; i++) {
        double k1 = h * f(option, x, y);
        double k2 = h * f(option, x + h/2.0, y + k1/2.0);
        double k3 = h * f(option, x + h/2.0, y + k2/2.0);
        double k4 = h * f(option, x + h,     y + k3);

        cout << left
             << setw(6)  << step
             << setw(12) << x
             << setw(14) << y
             << setw(14) << k1
             << setw(14) << k2
             << setw(14) << k3
             << setw(14) << k4
             << "\n";

        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        x = x + h;
        step++;
    }

    cout << "\nApproximate solution at x = " << x << " is y = " << y << "\n";
    return 0;
}

```

---

## RK4 Input
**Input (input.txt):**   
```
3
0 0.5
0.2 0.1
```
---

## Output
**Output (output.txt):** 
```
Runge-Kutta 4th Order (RK4)
Chosen function option = 3
Initial condition: x0 = 0.000000, y0 = 0.500000
Target: x = 0.200000, step size h = 0.100000

Step  x           y             k1            k2            k3            k4            
----------------------------------------------------------------------------------------
0     0.000000    0.500000      0.150000      0.157250      0.157613      0.164761      
1     0.100000    0.657414      0.164741      0.171729      0.172078      0.178949      

Approximate solution at x = 0.200000 is y = 0.829298
```
---
---

# Curve Fitting
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-red?style=for-the-badge)](./Curve%20Fitting/)
## 📐 Mathematical Foundation
Curve fitting is the process of finding a mathematical equation that best represents a given set of data points.
The most common approach is the **Least Squares Method**, which minimizes the sum of squares of the errors
between observed values and computed values obtained from the fitted curve.

You will find implementations for:
- Linear curve fitting using the **Least Squares Line**
- Polynomial curve fitting using the **Least Squares Polynomial**
- **Non-linear curve fitting** using suitable mathematical transformations

Each method includes:
- Theory and algorithm explanation
- C++ source code
- Sample input and output files


### 🔍 Overview

| Method | Curve Type | Main Idea | Notes |
|---|---|---|---|
| Least Squares Line | Linear | Minimize squared errors | Simple and widely used |
| Least Squares Polynomial | Quadratic | Higher-degree fitting | Better for curved data |
| Non-Linear Curve Fitting | Exponential / Power | Transform to linear form | Model-dependent accuracy |
---

# Least Squares Line
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Curve%20Fitting/Least%20Squares%20Line/)

## Least Squares Line Theory
### 📌 Mathematical Model
The straight line equation is assumed as:

y = a + b·x

where a and b are constants to be determined.


### 📐 Normal Equations
Using the least squares principle, the normal equations are:

sum(y) = n·a + b·sum(x)

sum(xy) = a·sum(x) + b·sum(x²)

Solving these equations gives the values of a and b.


### 🧾 Algorithm
1. Read the number of data points n
2. Read values of x and y
3. Compute sum(x), sum(y), sum(xy), sum(x²)
4. Form the normal equations
5. Solve for a and b
6. Display the fitted line equation
---

## Least Squares Line Code
```cpp

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

int main() {
    int n;
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin || !fout) {
        cerr << "Error opening file." << endl;
        return 1;
    }

    fin >> n;
    double x[100], y[100];
    double sumX=0, sumY=0, sumXY=0, sumX2=0;

    fout << fixed << setprecision(6);
    fout << "Number of points: " << n << endl;
    fout << "\nInput Points (x, y):" << endl;
    for(int i=0;i<n;i++){
        fin >> x[i] >> y[i];
        fout << "(" << x[i] << ", " << y[i] << ")" << endl;
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i]*y[i];
        sumX2 += x[i]*x[i];
    }

    fout << "\nSum of x: " << sumX << endl;
    fout << "Sum of y: " << sumY << endl;
    fout << "Sum of x*y: " << sumXY << endl;
    fout << "Sum of x^2: " << sumX2 << endl;

    double b = (n*sumXY - sumX*sumY) / (n*sumX2 - sumX*sumX);
    double a = (sumY - b*sumX) / n;

    fout << "\nEquation of best fit line (Least Squares):" << endl;
    fout << "y = " << a << " + " << b << "x" << endl;

    fin.close();
    fout.close();
    return 0;
}
```

---

## Least Squares Line Input
**Input (input.txt):**   
```
5
1 2
2 3
3 5
4 4
5 6
```
---

## Least Squares Line Output
**Output (output.txt):** 
```
Number of points: 5

Input Points (x, y):
(1.000000, 2.000000)
(2.000000, 3.000000)
(3.000000, 5.000000)
(4.000000, 4.000000)
(5.000000, 6.000000)

Sum of x: 15.000000
Sum of y: 20.000000
Sum of x*y: 69.000000
Sum of x^2: 55.000000

Equation of best fit line (Least Squares):
y = 1.300000 + 0.900000x
```
---

# Least Squares Polynomial
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Curve%20Fitting/Least%20Squares%20Polynomial/)

## Least Squares Polynomial Theory
### 📌 Mathematical Model
The quadratic polynomial is assumed as:

y = a + b·x + c·x²

where a, b, and c are constants.


### 📐 Normal Equations
The normal equations are:

sum(y) = n·a + b·sum(x) + c·sum(x²)

sum(xy) = a·sum(x) + b·sum(x²) + c·sum(x³)

sum(x²y) = a·sum(x²) + b·sum(x³) + c·sum(x⁴)

Solving these equations gives the polynomial coefficients.


### 🧾 Algorithm
1. Read the number of observations
2. Read the values of x and y
3. Compute required summations
4. Form the normal equations
5. Solve for a, b, and c
6. Display the fitted polynomial equation

---

## Least Squares Polynomial Code
```cpp

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

int main() {
    int n;
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin || !fout) {
        cerr << "Error opening file." << endl;
        return 1;
    }

    fin >> n;
    double x[100], y[100];
    fout << fixed << setprecision(6);
    fout << "Number of points: " << n << endl;
    fout << "\nInput Points (x, y):" << endl;
    for(int i=0;i<n;i++) {
        fin >> x[i] >> y[i];
        fout << "(" << x[i] << ", " << y[i] << ")" << endl;
    }

    double sumX=0,sumX2=0,sumX3=0,sumX4=0,sumY=0,sumXY=0,sumX2Y=0;
    for(int i=0;i<n;i++){
        sumX += x[i];
        sumX2 += x[i]*x[i];
        sumX3 += x[i]*x[i]*x[i];
        sumX4 += x[i]*x[i]*x[i]*x[i];
        sumY += y[i];
        sumXY += x[i]*y[i];
        sumX2Y += x[i]*x[i]*y[i];
    }

    fout << "\nSum of x: " << sumX << endl;
    fout << "Sum of x^2: " << sumX2 << endl;
    fout << "Sum of x^3: " << sumX3 << endl;
    fout << "Sum of x^4: " << sumX4 << endl;
    fout << "Sum of y: " << sumY << endl;
    fout << "Sum of x*y: " << sumXY << endl;
    fout << "Sum of x^2*y: " << sumX2Y << endl;

    // Solve normal equations for quadratic fit: y = a + b*x + c*x^2
    // [ n    sumX   sumX2 ] [a]   [ sumY   ]
    // [sumX sumX2  sumX3 ] [b] = [ sumXY  ]
    // [sumX2 sumX3 sumX4] [c]   [ sumX2Y ]
    double A[3][4] = {
        {double(n), sumX, sumX2, sumY},
        {sumX, sumX2, sumX3, sumXY},
        {sumX2, sumX3, sumX4, sumX2Y}
    };

    // Gaussian elimination
    for(int i=0;i<3;i++){
        // Partial pivoting
        int maxRow = i;
        for(int k=i+1;k<3;k++){
            if(abs(A[k][i]) > abs(A[maxRow][i])) maxRow = k;
        }
        for(int k=i;k<4;k++) swap(A[maxRow][k], A[i][k]);
        // Eliminate
        for(int k=i+1;k<3;k++){
            double f = A[k][i]/A[i][i];
            for(int j=i;j<4;j++)
                A[k][j] -= f*A[i][j];
        }
    }
    // Back substitution
    double coeff[3];
    for(int i=2;i>=0;i--){
        coeff[i] = A[i][3];
        for(int j=i+1;j<3;j++)
            coeff[i] -= A[i][j]*coeff[j];
        coeff[i] /= A[i][i];
    }

    fout << "\nEquation of best fit quadratic polynomial (Least Squares):" << endl;
    fout << "y = " << coeff[0] << " + " << coeff[1] << "x + " << coeff[2] << "x^2" << endl;

    fin.close();
    fout.close();
    return 0;
}
```

---

## Least Squares Polynomial Input
**Input (input.txt):**   
```
4
1 1
2 4
3 9
4 16
```
---

## Least Squares Polynomial Output
**Output (output.txt):** 
```
Number of points: 4

Input Points (x, y):
(1.000000, 1.000000)
(2.000000, 4.000000)
(3.000000, 9.000000)
(4.000000, 16.000000)

Sum of x: 10.000000
Sum of x^2: 30.000000
Sum of x^3: 100.000000
Sum of x^4: 354.000000
Sum of y: 30.000000
Sum of x*y: 100.000000
Sum of x^2*y: 354.000000

Equation of best fit quadratic polynomial (Least Squares):
y = 0.000000 + -0.000000x + 1.000000x^2
```
---

# Non Linear Curve Fitting
[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-purple?style=for-the-badge)](./Curve%20Fitting/Non%20Linear%20Curve%20Fitting/)

## Non Linear Curve Fitting Theory
### 📌 Common Models
Some widely used non-linear models are:

- **Exponential curve**
y = a·e^(b·x)

- **Power curve**
y = a·x^b


### 🔄 Linearization Technique
For the exponential model:

y = a·e^(b·x)

Taking logarithm on both sides:

ln(y) = ln(a) + b·x

This converts the equation into a linear form suitable for least squares fitting.

### 🧾 Algorithm
1. Read the given data points
2. Apply logarithmic transformation
3. Convert the equation into linear form
4. Apply least squares method
5. Compute constants
6. Convert back to original non-linear equation
7. Display the fitted curve
---

## Non Linear Curve Fitting Code
```cpp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

int main() {
    int n;
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin || !fout) {
        cerr << "Error opening file." << endl;
        return 1;
    }

    fin >> n;
    double x[100], y[100];
    fout << fixed << setprecision(6);
    fout << "Number of points: " << n << endl;
    fout << "\nInput Points (x, y):" << endl;
    for(int i=0;i<n;i++) {
        fin >> x[i] >> y[i];
        fout << "(" << x[i] << ", " << y[i] << ")" << endl;
    }

    // Fit y = a * exp(bx) using linearization: ln(y) = ln(a) + b*x
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    double lnY[100];
    for(int i=0;i<n;i++) {
        lnY[i] = log(y[i]);
        sumX += x[i];
        sumY += lnY[i];
        sumXY += x[i]*lnY[i];
        sumX2 += x[i]*x[i];
    }

    fout << "\nSum of x: " << sumX << endl;
    fout << "Sum of ln(y): " << sumY << endl;
    fout << "Sum of x*ln(y): " << sumXY << endl;
    fout << "Sum of x^2: " << sumX2 << endl;

    double b = (n*sumXY - sumX*sumY) / (n*sumX2 - sumX*sumX);
    double A = (sumY - b*sumX) / n;
    double a = exp(A);

    fout << "\nEquation of best fit (Non-linear, y = a*exp(bx)):" << endl;
    fout << "y = " << a << " * exp(" << b << "x)" << endl;

    fin.close();
    fout.close();
    return 0;
}

```

---

## Non Linear Curve Fitting Input
**Input (input.txt):**   
```
4
1 2.7
2 7.4
3 20.1
4 54.6
```
---

## Non Linear Curve Fitting Output
**Output (output.txt):** 
```
Number of points: 4

Input Points (x, y):
(1.000000, 2.700000)
(2.000000, 7.400000)
(3.000000, 20.100000)
(4.000000, 54.600000)

Sum of x: 10.000000
Sum of ln(y): 9.995485
Sum of x*ln(y): 29.998507
Sum of x^2: 30.000000

Equation of best fit (Non-linear, y = a*exp(bx)):
y = 0.993993 * exp(1.001959x)
```

---
---

<a id="authors"></a>
# 👨‍💻 Authors

## Project Lead
### Abir Hasan Arko
[![GitHub](https://img.shields.io/badge/GitHub-AbirHasanArko-181717?style=flat&logo=github)](https://github.com/AbirHasanArko)  
CSE, KUET  
Roll: 2207053

### Contributions:   
[![Solution of Linear Equations](https://img.shields.io/badge/📂-Solution%20of%20Linear%20Equations-blue?style=for-the-badge)](./Solution%20of%20Linear%20Equations/)  
[![Interpolation and Approximation](https://img.shields.io/badge/📂-Interpolation%20and%20Approximation-brown?style=for-the-badge)](./Interpolation%20and%20Approximation/)  
[![Numerical Differentiation](https://img.shields.io/badge/📂-Numerical%20Differentiation-purple?style=for-the-badge)](./Numerical%20Differentiation/)   
[![Curve Fitting](https://img.shields.io/badge/📂-Curve%20Fitting-green?style=for-the-badge)](./Curve%20Fitting/)  

---

## Principal contributor
### Ajoy Saha
[![GitHub](https://img.shields.io/badge/GitHub-ajoycodes-181717?style=flat&logo=github)](https://github.com/ajoycodes)  
CSE, KUET  
Roll: 2207037

### Contributions:   
[![Numerical Integration](https://img.shields.io/badge/📂-Numerical%20Integration-blue?style=for-the-badge)](./Numerical%20Integration/)  
[![Solution of Differential Equations](https://img.shields.io/badge/📂-Solution%20of%20Differential%20Equations-brown?style=for-the-badge)](./Solution%20of%20Differential%20Equations/)  
[![Curve Fitting](https://img.shields.io/badge/📂-Curve%20Fitting-purple?style=for-the-badge)](./Curve%20Fitting/)  

---

## Principal contributor
### Md. Shomik Shahriar
[![GitHub](https://img.shields.io/badge/GitHub-HapiGuy-181717?style=flat&logo=github)](https://github.com/Hapi-Guy)  
CSE, KUET  
Roll: 2207041

### Contributions:   
[![Solution of Non-linear Equations](https://img.shields.io/badge/📂-Solution%20of%20Non%20Linear%20Equations-blue?style=for-the-badge)](./Solution%20of%20Non-linear%20Equations/)  
[![Numerical Differentiation](https://img.shields.io/badge/📂-Numerical%20Differentiation-purple?style=for-the-badge)](./Numerical%20Differentiation/)  

---
