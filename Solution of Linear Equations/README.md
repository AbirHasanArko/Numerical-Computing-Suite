# Solution of Linear Equations

[![C++](https://img.shields.io/badge/Language-C++-00599C?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)
[![Numerical Methods](https://img.shields.io/badge/Topic-Numerical%20Methods-FF6B6B?style=for-the-badge)](https://github.com/AbirHasanArko/Numerical-Computing-Suite)

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Overview of Solution Methods](#-overview-of-solution-methods)
- [Mathematical Foundation](#-mathematical-foundation)
  - [What are Linear Equations?](#what-are-linear-equations)
  - [Matrix Representation](#matrix-representation)
  - [Types of Solutions](#types-of-solutions)
  - [Fundamental Concepts](#fundamental-concepts)
- [Solution Methods](#-solution-methods)
  - [1. Gauss Elimination Method](#1-gauss-elimination-method)
  - [2. Gauss-Jordan Elimination Method](#2-gauss-jordan-elimination-method)
  - [3. LU Decomposition Method](#3-lu-decomposition-method)
- [Method Comparison](#-method-comparison)
- [Applications](#-applications)
- [Implementation Structure](#-implementation-structure)
- [Getting Started](#-getting-started)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

This section of the **Numerical Computing Suite** focuses on solving **systems of linear equations** using various numerical methods. A system of linear equations is a collection of linear equations involving the same set of variables, and finding solutions to these systems is fundamental to many areas of mathematics, science, and engineering.

Linear systems appear in countless applications: 
- **Engineering**: Circuit analysis, structural analysis, control systems
- **Physics**: Equilibrium problems, quantum mechanics
- **Economics**: Input-output models, optimization
- **Computer Graphics**:  Transformations, rendering
- **Data Science**: Linear regression, machine learning algorithms

This collection provides three powerful methods for solving linear systems, each with its own advantages and use cases.  All implementations include detailed step-by-step visualization, solution type detection, and comprehensive error handling. 

---

## 🔍 Overview of Solution Methods

This repository implements three classical numerical methods:  

| Method | Best Used For | Complexity | Key Feature |
|--------|--------------|------------|-------------|
| **Gauss Elimination** | Single solution, general systems | O(n³) | Forward elimination + back substitution |
| **Gauss-Jordan Elimination** | Finding inverse, reduced form | O(n³) | Complete diagonal reduction |
| **LU Decomposition** | Multiple systems, repeated solving | O(n³) | Matrix factorization |

---

## 📐 Mathematical Foundation

### What are Linear Equations?

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

## 🛠️ Solution Methods

### 1. Gauss Elimination Method

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-blue?style=for-the-badge)](./Gauss%20Elimination%20Method/)

#### Theory

**Gauss Elimination** (named after Carl Friedrich Gauss, 1777-1855) is a systematic algorithm for solving systems of linear equations.  It's one of the most important and widely-used methods in numerical linear algebra, forming the foundation for many advanced techniques.

#### Historical Context

Though named after Gauss, the method was known to ancient Chinese mathematicians (documented in "The Nine Chapters on the Mathematical Art" around 150 BC). Gauss popularized it in the West and extended its applications to least squares problems and geodetic calculations.

#### Mathematical Foundation

The method is based on the principle that **elementary row operations** preserve the solution set of a linear system. These operations are:

1. **Row Swapping**: Interchange two rows (R₁ ↔ R₂)
2. **Row Scaling**: Multiply a row by a non-zero scalar (R₁ → kR₁)
3. **Row Addition**: Add a multiple of one row to another (R₂ → R₂ + kR₁)

**Key Insight**: These operations transform the system into an equivalent (same solution) but simpler form.

#### The Two-Phase Approach

##### **Phase 1: Forward Elimination**

**Goal**: Transform the augmented matrix [A|b] into **upper triangular form** (row echelon form).

**Process**:  For each column k from 1 to n-1:

1. **Pivot Selection** (with partial pivoting):
   - Find the row i ≥ k with the largest |aᵢₖ| in column k
   - Swap row k with row i
   - This element aₖₖ becomes the **pivot**

2. **Elimination Step**:
   - For each row i below the pivot (i = k+1 to n):
     - Compute multiplier:  `mᵢₖ = aᵢₖ / aₖₖ`
     - Update row i: `Row_i = Row_i - mᵢₖ × Row_k`
     - This makes aᵢₖ = 0 (eliminates the element below the pivot)

**Mathematical Representation**: 

Starting with:
```
┌                    ┐
│ a₁₁  a₁₂  a₁₃ | b₁ │
│ a₂₁  a₂₂  a₂₃ | b₂ │
│ a₃₁  a₃₂  a₃₃ | b₃ │
└                    ┘
```

After forward elimination:
```
┌                     ┐
│ u₁₁  u₁₂  u₁₃ | b₁' │
│  0   u₂₂  u₂₃ | b₂' │
│  0    0   u₃₃ | b₃' │
└                     ┘
```

**Example Walkthrough**: 

Initial system:
```
2x + y - z = 8       ┌              ┐
-3x - y + 2z = -11   │  2  1 -1 | 8 │
-2x + y + 2z = -3    │ -3 -1  2 |-11│
                     │ -2  1  2 |-3 │
                     └              ┘
```

Step 1: Eliminate x from rows 2 and 3
- m₂₁ = -3/2 = -1.5
- m₃₁ = -2/2 = -1. 0

```
┌                    ┐
│  2   1  -1 |  8    │
│  0  0. 5 0.5 |  1  │  (R₂ - m₂₁×R₁)
│  0   2   1  |  5   │  (R₃ - m₃₁×R₁)
└                    ┘
```

Step 2: Eliminate y from row 3
- m₃₂ = 2/0.5 = 4

```
┌                    ┐
│  2  1   -1 |  8    │
│  0  0.5 0.5|  1    │
│  0  0   -1 |  1    │  (R₃ - m₃₂×R₂)
└                    ┘
```

##### **Phase 2: Back Substitution**

**Goal**:  Solve for variables starting from the last equation and working backwards.

**Process**: 

1. **Last Variable**: From the last row, solve directly
   ```
   xₙ = b'ₙ / uₙₙ
   ```

2. **Remaining Variables**: For i = n-1, n-2, ..., 1:
   ```
   xᵢ = (b'ᵢ - Σⱼ₌ᵢ₊₁ⁿ uᵢⱼ × xⱼ) / uᵢᵢ
   ```

**Example (continued)**:

From the upper triangular form:
```
2x + y - z = 8
0.5y + 0.5z = 1
-z = 1
```

Solve backwards:
1. z = -1 (from equation 3)
2. y = (1 - 0.5(-1))/0.5 = 3 (substitute z into equation 2)
3. x = (8 - 3 - (-1))/2 = 2 (substitute y and z into equation 1)

**Solution**: x = 2, y = 3, z = -1 ✓

#### Partial Pivoting Strategy

**Why Pivoting?**

Consider this system without pivoting:
```
0.0001x + y = 1
x + y = 2
```

If we eliminate using the small pivot 0.0001:
- Multiplier = 1/0.0001 = 10,000
- This magnifies round-off errors dramatically! 

**Partial Pivoting Solution**: 
- Always choose the largest available pivot in the current column
- Swap rows to bring it to the pivot position
- This minimizes error propagation

**Algorithm**: 
```
For each column k: 
    Find i ≥ k such that |aᵢₖ| is maximum
    Swap row k with row i
    Proceed with elimination
```

**Benefits**:
- ✅ Reduces round-off error accumulation
- ✅ Avoids division by small numbers
- ✅ Improves numerical stability
- ✅ Often prevents division by zero

**Trade-off**: Adds O(n²) comparisons, but the stability gain is worth it.

#### Solution Type Detection

After forward elimination, analyze the resulting upper triangular matrix:

**Step 1: Calculate Rank**
- Count non-zero rows (rows where at least one coefficient ≠ 0)
- This gives rank(A) after elimination

**Step 2: Check Consistency**
- If any row has form [0 0 ... 0 | c] where c ≠ 0:
  - This means 0 = c, which is impossible
  - System is **inconsistent** → **No Solution**

**Step 3: Compare Ranks**
- If rank(A) = rank([A|b]) = n:
  - System is consistent and determined → **Unique Solution**
- If rank(A) = rank([A|b]) < n:
  - System is consistent but underdetermined → **Infinite Solutions**

**Visual Decision Tree**:
```
After Forward Elimination
        |
        ├─→ Found [0 0 ... 0 | c≠0]?  → YES → No Solution ❌
        |
        └─→ NO → rank = n? ─→ YES → Unique Solution ✓
                      |
                      └─→ NO → Infinite Solutions ∞
```

#### Complexity Analysis

**Time Complexity**: 

1. **Forward Elimination**:
   - Outer loop: k = 1 to n-1 (n-1 iterations)
   - For each k: 
     - Pivoting: O(n) comparisons
     - For each row i > k: O(n) operations
   - Total:  Σₖ₌₁ⁿ⁻¹ (n-k)² ≈ n³/3 operations

2. **Back Substitution**:
   - For variable i: need to compute i-1 multiplications and additions
   - Total:  Σᵢ₌₁ⁿ i = n(n+1)/2 ≈ n²/2 operations

**Overall**:  O(n³) dominated by forward elimination

**Space Complexity**:  O(n²) for storing the augmented matrix

**Exact Operation Counts** (for n×n system):
- **Multiplications/Divisions**: (2n³ + 3n² - 5n)/6 ≈ n³/3
- **Additions/Subtractions**: (n³ - n)/3 ≈ n³/3
- **Comparisons** (with pivoting): n(n-1)/2 ≈ n²/2

#### Numerical Considerations

**Round-off Error Sources**:
1. **Subtractive Cancellation**: When subtracting nearly equal numbers
2. **Small Pivots**: Division by small numbers amplifies errors
3. **Error Propagation**: Errors in early steps affect later computations

**Mitigation Strategies**:
- Use **partial pivoting** (or complete pivoting for extreme cases)
- Use **double precision** floating-point arithmetic
- **Equilibration**:  Scale rows/columns to have similar magnitudes
- **Iterative refinement**: Use solution to improve accuracy

#### Advantages & Limitations

**Advantages**:  
✅ Straightforward and intuitive algorithm  
✅ Efficient for single right-hand side  
✅ Foundation for many other methods  
✅ Works for any non-singular square system  
✅ Easy to implement and debug  
✅ Numerical stability with pivoting  

**Limitations**:  
❌ Not efficient for multiple right-hand sides (must repeat for each)  
❌ Requires O(n²) storage  
❌ Sensitive to round-off errors without pivoting  
❌ Inefficient for sparse matrices (many zeros)  
❌ Cannot exploit special matrix structures  

---

### 2. Gauss-Jordan Elimination Method

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-green?style=for-the-badge)](./Gauss-Jordan%20Elimination%20Method/)

#### Theory

**Gauss-Jordan Elimination** is an extension of Gaussian elimination that reduces the augmented matrix to **reduced row echelon form (RREF)** rather than just row echelon form.  Named after Carl Friedrich Gauss and Wilhelm Jordan (1842-1899), this method produces the solution directly without requiring back substitution.

#### Historical Context

Wilhelm Jordan, a German geodesist, popularized this variant in his 1888 handbook on geodesy. While Gauss elimination stops at triangular form, Jordan's method continues the elimination process to achieve complete diagonal reduction, making it particularly useful for matrix inversion.

#### Mathematical Foundation

The method extends the three elementary row operations to achieve a more complete reduction:

**Goal**: Transform [A|b] into [I|x], where: 
- **I** is the identity matrix
- **x** is the solution vector

**Key Difference from Gauss Elimination**: 
- Gauss:  Eliminate only **below** each pivot
- Gauss-Jordan: Eliminate **both above and below** each pivot

#### The Complete Reduction Process

##### **Phase 1: Forward Pass**

Similar to Gauss elimination, but with an additional step: 

For each column k from 1 to n: 

1. **Partial Pivoting**:
   - Find row i ≥ k with maximum |aᵢₖ|
   - Swap row k with row i

2. **Pivot Normalization**:
   - Divide entire row k by aₖₖ to make pivot = 1
   - `Row_k = Row_k / aₖₖ`

3. **Complete Column Elimination**:
   - For **ALL** rows i ≠ k (not just i > k):
     - Compute multiplier: `mᵢₖ = aᵢₖ`
     - Update:  `Row_i = Row_i - mᵢₖ × Row_k`
   - This makes all elements in column k (except the pivot) equal to zero

**Mathematical Representation**: 

Starting augmented matrix:
```
┌                    ┐
│ a₁₁  a₁₂  a₁₃ | b₁ │
│ a₂₁  a₂₂  a₂₃ | b₂ │
│ a₃₁  a₃₂  a₃₃ | b₃ │
└                    ┘
```

After Gauss-Jordan elimination:
```
┌                      ┐
│  1    0    0  | x₁   │
│  0    1    0  | x₂   │
│  0    0    1  | x₃   │
└                      ┘
```

The solution is **directly visible**:  x₁, x₂, x₃ are in the last column! 

#### Detailed Example Walkthrough

**Initial System**: 
```
2x + y - z = 8
-3x - y + 2z = -11
-2x + y + 2z = -3
```

**Augmented Matrix**:
```
┌              ┐
│  2  1 -1 | 8 │
│ -3 -1  2 |-11│
│ -2  1  2 |-3 │
└              ┘
```

**Step 1: Process column 1**

a) Pivot is at (1,1): pivot = 2 (or swap if needed)

b) Normalize row 1 (divide by 2):
```
┌                 ┐
│  1  0.5 -0.5|4  │  (R₁/2)
│ -3  -1    2 |-11│
│ -2   1    2 |-3 │
└                 ┘
```

c) Eliminate column 1 in rows 2 and 3:
```
┌                   ┐
│  1  0.5 -0.5  | 4 │
│  0  0. 5  0.5 | 1 │  (R₂ + 3R₁)
│  0  2    1    | 5 │  (R₃ + 2R₁)
└                   ┘
```

**Step 2: Process column 2**

a) Pivot is at (2,2): pivot = 0. 5

b) Normalize row 2 (divide by 0.5):
```
┌                   ┐
│  1  0.5 -0.5 | 4  │
│  0  1    1   | 2  │  (R₂/0.5)
│  0  2    1   | 5  │
└                   ┘
```

c) Eliminate column 2 in rows 1 and 3:
```
┌                   ┐
│  1  0   -1   | 3  │  (R₁ - 0.5R₂)
│  0  1    1   | 2  │
│  0  0   -1   | 1  │  (R₃ - 2R₂)
└                   ┘
```

**Step 3: Process column 3**

a) Pivot is at (3,3): pivot = -1

b) Normalize row 3 (divide by -1):
```
┌                   ┐
│  1  0  -1   | 3   │
│  0  1   1   | 2   │
│  0  0   1   |-1   │  (R₃/-1)
└                   ┘
```

c) Eliminate column 3 in rows 1 and 2:
```
┌                   ┐
│  1  0   0   | 2   │  (R₁ + R₃)
│  0  1   0   | 3   │  (R₂ - R₃)
│  0  0   1   |-1   │
└                   ┘
```

**Final Result**:  [I|x] form achieved!

**Solution reads directly**: x = 2, y = 3, z = -1 ✓

#### Key Advantages of RREF

**1. Direct Solution Reading**:
- No back substitution needed
- Solution is immediately visible in the last column
- Less prone to calculation errors in the final step

**2. Matrix Inversion**:
To find A⁻¹, augment A with identity:   [A|I]
Apply Gauss-Jordan:  Result is [I|A⁻¹]

**Example**: Find inverse of 2×2 matrix
```
[2 1|1 0]     [1 0. 5|0.5   0]     [1 0|0.5  -0.5]
[1 3|0 1]  →  [1   3|  0   1]  →  [0 1|-0.2  0.4]
```
Inverse is: 
```
A⁻¹ = [ 0.5  -0.5]
      [-0.2   0.4]
```

**3. Rank Determination**:
- The number of non-zero rows in RREF equals the rank
- More explicit than row echelon form

**4. Solution Space Visualization**:
- For infinite solutions, free variables are immediately identifiable
- Parametric solution form is easier to construct

#### Comparison with Standard Gauss Elimination

| Aspect | Gauss Elimination | Gauss-Jordan |
|--------|------------------|--------------|
| **Final Form** | Upper triangular | Diagonal (Identity) |
| **Pivots** | Can be any non-zero | Always 1 |
| **Elimination** | Below pivot only | Above & below pivot |
| **Back Substitution** | Required | Not needed |
| **Operations** | ~n³/3 | ~n³/2 |
| **Best For** | Single solution | Matrix inverse |
| **Intuitive** | Moderate | Very intuitive |

#### Algorithm Complexity

**Time Complexity**: 

1. **Forward Pass with Elimination**:
   - For each column k (n iterations):
     - Normalization: O(n) operations
     - Elimination in (n-1) rows: O(n²) operations
   - Total: n × n² = n³ operations

**Exact Count**: ~n³/2 operations (about 50% more than Gauss)

**Space Complexity**: O(n²) for the augmented matrix

**Why More Operations?**
- Gauss:  Eliminates only below pivot → triangular work
- Gauss-Jordan:  Eliminates above and below → rectangular work

**Operation Breakdown**:
- **Multiplications/Divisions**: ~n³/2
- **Additions/Subtractions**: ~n³/2
- **Total**: ~n³ (approximately 1.5 times Gauss elimination)

#### Numerical Stability

**Stability Considerations**: 

1. **Pivoting is Essential**:
   - Without pivoting, method can be highly unstable
   - Partial pivoting is standard practice
   - Complete pivoting rarely needed but possible

2. **Error Propagation**:
   - More elimination steps → more opportunities for error
   - Errors in early steps affect more later steps
   - Partial pivoting mitigates this significantly

3. **Conditioning**:
   - Method's accuracy depends on matrix condition number κ(A)
   - Well-conditioned:  κ(A) ≈ 1 → reliable results
   - Ill-conditioned: κ(A) >> 1 → potential accuracy loss

**Condition Number**:
```
κ(A) = ||A|| × ||A⁻¹||
```

- κ(A) = 1: Perfect conditioning (orthogonal matrices)
- κ(A) < 10³: Well-conditioned
- κ(A) > 10⁶: Ill-conditioned (problematic)

#### Special Applications

**1. Finding Matrix Inverse**:
```
Start:  [A|I]
End:    [I|A⁻¹]
```

**2. Solving Multiple Systems**: 
Augment with multiple b vectors:  [A|b₁ b₂ ... bₖ]
Result: [I|x₁ x₂ ...  xₖ]

**3. Rank Computation**:
The number of pivots (non-zero rows in RREF) = rank(A)

**4. Basis Finding**: 
Pivot columns in original matrix form a basis for column space

**5. Null Space**:
Free variables in RREF reveal null space basis vectors

#### When to Use Gauss-Jordan

**Optimal Scenarios**:  
✅ Finding matrix inverse  
✅ Solving Ax = b for multiple b vectors simultaneously  
✅ Pedagogical purposes (teaching linear algebra)  
✅ Explicit RREF needed for analysis  
✅ When back substitution code is error-prone  
✅ Small to medium-sized dense matrices  

**Avoid When**:  
❌ Solving single system (Gauss elimination is faster)  
❌ Very large matrices (extra operations matter)  
❌ Sparse matrices (destroys sparsity pattern)  
❌ When only solution needed (not RREF)  
❌ Numerical stability is critical (LU better)  

#### Advantages & Limitations

**Advantages**:  
✅ No back substitution required  
✅ Solution directly readable  
✅ Excellent for matrix inversion  
✅ Very intuitive and teachable  
✅ Handles multiple right-hand sides efficiently  
✅ Clear geometric interpretation  

**Limitations**:   
❌ More operations than standard Gauss (~50% more)  
❌ Not optimal for single system  
❌ More floating-point operations → more round-off error  
❌ Destroys matrix sparsity  
❌ Inefficient for very large systems  

---

### 3. LU Decomposition Method

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-orange?style=for-the-badge)](./LU%20Decomposition/)

#### Theory

**LU Decomposition** (also called LU Factorization) is a matrix factorization method that decomposes a square matrix **A** into the product of a **Lower triangular matrix (L)** and an **Upper triangular matrix (U)**: 

```
A = L × U
```

This powerful technique transforms the problem of solving Ax = b into solving two simpler triangular systems.  It's one of the most important factorizations in numerical linear algebra and forms the basis for many advanced algorithms.

#### Historical Context

The systematic approach to LU decomposition was developed in the 1940s by several mathematicians including mathematician Alan Turing.  The method became practical with the advent of computers, as it allows efficient solutions of linear systems, especially when the matrix A remains fixed but the right-hand side b changes multiple times.

#### Mathematical Foundation

**Core Principle**: Instead of solving Ax = b directly, factor A first: 

```
A = LU
Ax = b  becomes  LUx = b
```

**Two-Step Solution**:
1. Solve **Ly = b** for y (forward substitution)
2. Solve **Ux = y** for x (back substitution)

**Why This Works**:
- Triangular systems are easy to solve:  O(n²) time
- Factorization is done once: O(n³) time
- For k different b vectors:  O(n³ + kn²) instead of O(kn³)
- **Massive savings** when k > 1! 

#### Types of LU Decomposition

There are several variants of LU decomposition based on how L and U are defined:

##### **1. Doolittle's Method** (Used in our implementation)

**Definition**:
- L has **1's on the diagonal**
- U has **computed values on the diagonal**

**Form**:
```
A = L × U

┌             ┐   ┌           ┐   ┌             ┐
│ a₁₁ a₁₂ a₁₃ │   │ 1   0   0 │   │ u₁₁ u₁₂ u₁₃ │
│ a₂₁ a₂₂ a₂₃ │ = │ l₂₁ 1   0 │ × │ 0   u₂₂ u₂₃ │
│ a₃₁ a₃₂ a₃₃ │   │ l₃₁ l₃₂ 1 │   │ 0    0  u₃₃ │
└             ┘   └           ┘   └             ┘
```

**Characteristics**:
- Most commonly used
- Natural extension of Gaussian elimination
- U matrix is exactly what you get from forward elimination

##### **2. Crout's Method**

**Definition**:
- L has **computed values on the diagonal**
- U has **1's on the diagonal**

**Form**:
```
┌            ┐   ┌             ┐
│ l₁₁  0  0  │   │  1  u₁₂ u₁₃ │
│ l₂₁ l₂₂ 0  │ × │  0   1  u₂₃ │
│ l₃₁ l₃₂ l₃₃│   │  0   0   1  │
└            ┘   └             ┘
```

##### **3. Cholesky Decomposition** (For symmetric positive-definite matrices)

**Definition**: A = L × Lᵀ
- Special case where U = Lᵀ
- Only works for symmetric positive-definite matrices
- Requires ~n³/6 operations (half of standard LU)

#### Detailed Algorithm:  Doolittle's Method

##### **Step 1: Decomposition Process**

The decomposition proceeds column by column and row by row:

**For each column k from 1 to n:**

**Part A:  Compute U (Upper Triangular)**

For each element uᵢₖ where i ≤ k: 
```
uᵢₖ = aᵢₖ - Σⱼ₌₁ⁱ⁻¹ lᵢⱼ × uⱼₖ
```

**Part B: Compute L (Lower Triangular)**

For each element lₖᵢ where i > k:
```
lₖᵢ = (aₖᵢ - Σⱼ₌₁ᵏ⁻¹ lᵢⱼ × uⱼₖ) / uₖₖ
```

**Diagonal of L**: lₖₖ = 1 (by definition in Doolittle's method)

##### **Complete Algorithm** (Doolittle's Method)

```
Initialize L and U as n×n zero matrices
Set diagonal of L to 1

For i = 0 to n-1:
    
    // Compute U[i][k] for k = i to n-1
    For k = i to n-1:
        sum = 0
        For j = 0 to i-1:
            sum += L[i][j] × U[j][k]
        U[i][k] = A[i][k] - sum
    
    // Check for zero pivot
    If |U[i][i]| < ε:
        Matrix is singular (decomposition fails)
        Return error
    
    // Compute L[k][i] for k = i+1 to n-1
    For k = i+1 to n-1:
        sum = 0
        For j = 0 to i-1:
            sum += L[k][j] × U[j][i]
        L[k][i] = (A[k][i] - sum) / U[i][i]
```

#### Detailed Example:  3×3 Matrix

**Given System**:
```
2x + y - z = 8
-3x - y + 2z = -11
-2x + y + 2z = -3
```

**Matrix A**:
```
A = ┌           ┐
    │  2  1 -1  │
    │ -3 -1  2  │
    │ -2  1  2  │
    └           ┘
```

**Step-by-Step Decomposition**: 

**Iteration 1** (i = 0):

Compute U[0][k] for k = 0, 1, 2:
```
U[0][0] = A[0][0] = 2
U[0][1] = A[0][1] = 1
U[0][2] = A[0][2] = -1
```

Compute L[k][0] for k = 1, 2:
```
L[1][0] = A[1][0] / U[0][0] = -3 / 2 = -1.5
L[2][0] = A[2][0] / U[0][0] = -2 / 2 = -1.0
```

**Current State**:
```
L = ┌           ┐      U = ┌          ┐
    │  1  0  0  │          │ 2  1 -1  │
    │-1.5 1  0  │          │ 0  0  0  │
    │-1.0 0  1  │          │ 0  0  0  │
    └           ┘          └          ┘
```

**Iteration 2** (i = 1):

Compute U[1][k] for k = 1, 2:
```
U[1][1] = A[1][1] - L[1][0]×U[0][1]
        = -1 - (-1.5)×1
        = -1 + 1.5 = 0.5

U[1][2] = A[1][2] - L[1][0]×U[0][2]
        = 2 - (-1.5)×(-1)
        = 2 - 1.5 = 0.5
```

Compute L[2][1]: 
```
L[2][1] = (A[2][1] - L[2][0]×U[0][1]) / U[1][1]
        = (1 - (-1.0)×1) / 0.5
        = (1 + 1) / 0.5 = 4. 0
```

**Current State**:
```
L = ┌           ┐      U = ┌            ┐
    │  1   0  0 │          │ 2   1  -1  │
    │-1.5  1  0 │          │ 0  0.5 0.5 │
    │-1.0  4  1 │          │ 0   0   0  │
    └           ┘          └            ┘
```

**Iteration 3** (i = 2):

Compute U[2][2]: 
```
U[2][2] = A[2][2] - L[2][0]×U[0][2] - L[2][1]×U[1][2]
        = 2 - (-1.0)×(-1) - 4.0×0.5
        = 2 - 1 - 2 = -1.0
```

**Final Decomposition**:
```
L = ┌           ┐      U = ┌            ┐
    │  1   0  0 │          │ 2   1  -1  │
    │-1.5  1  0 │          │ 0  0.5 0.5 │
    │-1.0  4  1 │          │ 0   0 -1.0 │
    └           ┘          └            ┘
```

**Verification**:  Compute L × U
```
(L × U)[0][0] = 1×2 + 0×0 + 0×0 = 2 ✓
(L × U)[1][0] = -1.5×2 + 1×0 + 0×0 = -3 ✓
(L × U)[1][1] = -1.5×1 + 1×0.5 + 0×0 = -1 ✓
...  (all entries match A)
```

##### **Step 2: Forward Substitution (Solve Ly = b)**

Given:  Ly = b, find y

**Process**:  Solve from top to bottom
```
For i = 0 to n-1:
    y[i] = b[i]
    For j = 0 to i-1:
        y[i] -= L[i][j] × y[j]
    // Note: No division needed since L[i][i] = 1
```

**Example** (continuing from above, b = [8, -11, -3]):

```
y[0] = 8

y[1] = -11 - L[1][0]×y[0]
     = -11 - (-1.5)×8
     = -11 + 12 = 1

y[2] = -3 - L[2][0]×y[0] - L[2][1]×y[1]
     = -3 - (-1.0)×8 - 4.0×1
     = -3 + 8 - 4 = 1
```

**Result**: y = [8, 1, 1]

##### **Step 3: Back Substitution (Solve Ux = y)**

Given: Ux = y, find x

**Process**: Solve from bottom to top
```
For i = n-1 down to 0:
    x[i] = y[i]
    For j = i+1 to n-1:
        x[i] -= U[i][j] × x[j]
    x[i] /= U[i][i]
```

**Example** (y = [8, 1, 1]):

```
x[2] = 1 / U[2][2]
     = 1 / (-1.0) = -1

x[1] = (1 - U[1][2]×x[2]) / U[1][1]
     = (1 - 0.5×(-1)) / 0.5
     = (1 + 0.5) / 0.5 = 3

x[0] = (8 - U[0][1]×x[1] - U[0][2]×x[2]) / U[0][0]
     = (8 - 1×3 - (-1)×(-1)) / 2
     = (8 - 3 - 1) / 2 = 2
```

**Final Solution**: x = [2, 3, -1] ✓

#### The Power of LU:  Multiple Right-Hand Sides

**Scenario**: Solve Ax = b for different b vectors

**Without LU** (using Gauss elimination k times):
- Cost: k × (n³/3) ≈ kn³/3

**With LU**:
1.  Decompose A once: n³/3 operations
2. For each b:  Forward + Back substitution:  2n² operations
- Total cost: n³/3 + k(2n²) ≈ n³/3 + 2kn²

**Comparison** (for n = 1000):
- k = 1: LU ≈ same as Gauss
- k = 10: LU ≈ 90% faster! 
- k = 100: LU ≈ 99% faster!

**Example Application**: Circuit analysis with varying input voltages
- Circuit topology (matrix A) stays the same
- Input voltages (vector b) change
- Decompose once, solve many times efficiently! 

#### Connection to Gaussian Elimination

**Key Insight**: LU decomposition is essentially Gaussian elimination in disguise! 

**The Connection**:
- U is exactly the upper triangular matrix from forward elimination
- L encodes all the multipliers used during elimination

**Gauss Elimination**:
```
m₂₁ = a₂₁/a₁₁
Row2 = Row2 - m₂₁×Row1
```

**LU Decomposition**:
```
L[2][1] = m₂₁  (stores the multiplier)
U is the result after elimination
```

**Why Separate Them?**:
- LU explicitly stores the factorization
- Can reuse for different b vectors
- Enables additional operations (determinants, inverses)

#### Determinant Calculation

One beautiful property of LU decomposition: 

```
det(A) = det(L) × det(U)
       = 1 × det(U)           (since det(L) = 1 in Doolittle's)
       = ∏ᵢ₌₁ⁿ U[i][i]        (product of U's diagonal)
```

**Example**:
```
U = ┌            ┐
    │ 2   1  -1  │
    │ 0  0.5 0.5 │
    │ 0   0  -1  │
    └            ┘

det(A) = 2 × 0.5 × (-1) = -1
```

**Advantages**:
- Computing determinant is now O(n) after decomposition! 
- Without LU:  determinant computation is O(n³)

#### Pivoting in LU Decomposition

**The Problem**: A might need row swaps for stability

**Solution**: PA = LU (LU with Partial Pivoting)
- P is a **permutation matrix** (encodes row swaps)
- Solve:  PAx = Pb becomes LUx = Pb

**Algorithm Modification**:
```
At each step i:
    Find row k ≥ i with maximum |A[k][i]|
    Swap rows i and k in A
    Update permutation matrix P
    Continue with standard LU
```

**Permutation Matrix Example**:
```
P = ┌       ┐    Represents:  swap rows 1 and 2
    │ 0 1 0 │
    │ 1 0 0 │
    │ 0 0 1 │
    └       ┘
```

#### Complexity Analysis

**Time Complexity**: 

1. **LU Decomposition** (one-time cost):
   - Exact operations: (2n³ - 3n² + n)/6 ≈ n³/3
   - Dominant term: O(n³)

2. **Forward Substitution** (per solve):
   - Operations: n² - n ≈ n²
   - Complexity: O(n²)

3. **Back Substitution** (per solve):
   - Operations: n² ≈ n²
   - Complexity: O(n²)

**Total for k Right-Hand Sides**:
- Decomposition: ~n³/3
- k Solutions: ~2kn²
- **Overall**:  O(n³ + kn²)

**Comparison with Gauss Elimination** (k systems):
- Gauss: O(kn³/3)
- LU: O(n³/3 + 2kn²)
- **Crossover**: LU wins when k ≥ 2

**Space Complexity**:
- Store L, U, b, x: O(2n² + 2n) = O(n²)
- Can overwrite A with L and U: O(n²)

#### Numerical Stability and Conditioning

**Stability Factors**: 

1. **Growth Factor**:
   - Measures how large elements can become during decomposition
   - Without pivoting: can grow exponentially
   - With partial pivoting: usually bounded by 2ⁿ⁻¹ (rarely achieved)

2. **Condition Number Sensitivity**:
   - For ill-conditioned matrices (κ(A) >> 1):
     - Solution accuracy degrades
     - Pivoting essential
     - Consider iterative refinement

3. **Partial Pivoting Benefits**:
   - Usually sufficient for stability
   - Keeps maximum element in U ≤ 2×(max element in A)
   - Industry standard approach

**Error Bound**:
```
||x_computed - x_exact|| / ||x_exact|| ≈ κ(A) × machine_epsilon
```

#### Special Cases and Variants

**1. Tridiagonal Matrices**:
- Special LU algorithm:  O(n) time instead of O(n³)
- Common in differential equations

**2. Band Matrices**:
- Bandwidth b:  only O(nb²) operations
- Preserves band structure in L and U

**3. Sparse Matrices**:
- Fill-in problem: zeros can become non-zero
- Reordering strategies minimize fill-in
- Specialized sparse LU algorithms

**4. Positive Definite Matrices**: 
- Use Cholesky:  A = LLᵀ
- Half the operations of LU
- More stable

#### When to Use LU Decomposition

**Optimal Scenarios**:  
✅ **Solving multiple systems** with same A (k ≥ 2)  
✅ **Computing determinants** efficiently  
✅ **Matrix inversion** (solve Ax = eᵢ for each unit vector)  
✅ **Large systems** where reusability matters  
✅ **Numerical libraries** (LAPACK, BLAS implementations)  
✅ **Foundation for advanced methods** (iterative refinement)  

**Avoid When**:  
❌ Single system, small n (Gauss elimination simpler)  
❌ Sparse matrices without reordering  
❌ Ill-conditioned systems (use QR or SVD)  
❌ Memory constrained (stores L and U)  

#### Advantages & Limitations

**Advantages**:  
✅ **Extreme efficiency** for multiple right-hand sides  
✅ **One decomposition**, many solutions  
✅ **Easy determinant** calculation:  O(n) after decomposition  
✅ **Matrix inversion** natural application  
✅ **Foundation** for many advanced algorithms  
✅ **Numerical stability** with pivoting  
✅ **Explicit factorization** useful for analysis  
✅ **Parallelizable** to some extent  

**Limitations**:  
❌ **Higher memory** requirement (store L and U)  
❌ **Not optimal** for single solution  
❌ **Requires pivoting** for stability  
❌ **Fill-in problem** for sparse matrices  
❌ **Decomposition fails** if matrix is singular  
❌ **Not best** for ill-conditioned systems  
❌ **More complex** to implement correctly  

---

## 📊 Method Comparison

### When to Use Each Method?  

| Scenario | Recommended Method | Reason |
|----------|-------------------|--------|
| Solving one system once | Gauss Elimination | Simplest, least memory, adequate performance |
| Need solution without back substitution | Gauss-Jordan | Direct reading from RREF |
| Multiple systems, same A | **LU Decomposition** | Decompose once, massive savings |
| Finding matrix inverse | Gauss-Jordan or LU | Both efficient, GJ more straightforward |
| Numerical stability critical | LU with pivoting | Best error control |
| Teaching linear algebra | Gauss-Jordan | Most intuitive, clear steps |
| Large sparse matrices | Specialized methods | Standard methods destroy sparsity |
| Ill-conditioned systems | QR or SVD | Better numerical properties |

### Computational Cost Comparison

For solving `Ax = b` where A is n×n:

| Method | Single Solution | k Different b vectors | Memory |
|--------|----------------|---------------------|---------|
| **Gauss Elimination** | ~n³/3 + n² | ~k(n³/3) | O(n²) |
| **Gauss-Jordan** | ~n³/2 + n² | ~k(n³/2) | O(n²) |
| **LU Decomposition** | ~n³/3 + 2n² | ~n³/3 + 2kn² | O(2n²) |

**Concrete Example** (n = 1000):

| k | Gauss Ops | LU Ops | LU Speedup |
|---|-----------|---------|------------|
| 1 | 3.3×10⁸ | 3.3×10⁸ | ~1× |
| 10 | 3.3×10⁹ | 3.5×10⁸ | ~9. 4× |
| 100 | 3.3×10¹⁰ | 5.3×10⁸ | ~62× |

**Winner Analysis**:
- **k = 1**:  Gauss Elimination (slightly simpler)
- **k = 2-5**: LU starts winning
- **k > 5**: LU dramatically better

### Operation Count Details

**Forward Elimination/Decomposition**:
```
Gauss:         n³/3 - n²/2 + n/6
Gauss-Jordan:  n³/2
LU Doolittle:  n³/3
```

**Back/Forward Substitution**:
```
Gauss:        n²/2
Gauss-Jordan: 0 (included in elimination)
LU (both):    n²
```

### Stability Ranking

From most to least stable (with proper pivoting):

1. **LU with Partial Pivoting** ⭐⭐⭐⭐⭐
   - Industry standard
   - Excellent stability-performance balance

2. **Gauss with Partial Pivoting** ⭐⭐⭐⭐
   - Very stable
   - Slight edge to LU for multiple solves

3. **Gauss-Jordan with Pivoting** ⭐⭐⭐
   - More operations → more rounding errors
   - Still acceptable for most applications

**Note**: All methods MUST use pivoting for reliable results!

---

## 🎯 Applications

### Direct Applications

#### 1. **Circuit Analysis** ⚡
**Problem**: Find voltages and currents in electrical networks

**Approach**:  Apply Kirchhoff's laws
- Current law:  Σ currents = 0 at each node
- Voltage law:  Σ voltages = 0 around each loop
- Results in system of linear equations

**Why LU? **:  Circuit topology stays constant, but inputs vary
- Decompose circuit matrix once
- Solve for different voltage sources quickly

#### 2. **Structural Engineering** 🏗️
**Problem**: Analyze forces in trusses, beams, buildings

**Approach**: Equilibrium equations
- Force balance at each node
- Moment balance
- Linear system:  [Stiffness Matrix] × [Displacements] = [Forces]

**Why Gauss?**: Often single load case analysis

#### 3. **Chemical Engineering** ⚗️
**Problem**:  Material balance in reactors, distillation columns

**Approach**:  Conservation laws
- Mass balance
- Energy balance
- Component balance

**Why LU?**: Same process, different feed compositions

#### 4. **Economics** 💰
**Problem**:  Input-output models (Leontief models)

**Approach**: (I - A)x = d
- A: input-output matrix
- x: production levels
- d: final demand

**Why Gauss-Jordan?**: Often need to analyze matrix properties

### Computational Applications

#### 5. **Computer Graphics** 🎨
- **3D Transformations**: Solving for transformation matrices
- **Camera Calibration**: From known points to camera parameters
- **Lighting Calculations**:  Radiosity methods
- **Curve/Surface Fitting**: Interpolation control points

#### 6. **Machine Learning** 🤖
- **Linear Regression**: Normal equations Xᵀ Xx = Xᵀy
- **Principal Component Analysis**: Eigenvalue problems
- **Neural Network Training**: Weight updates
- **Support Vector Machines**: Quadratic programming subproblems

#### 7. **Physics Simulations** ⚛️
- **Quantum Mechanics**: Schrödinger equation discretization
- **Fluid Dynamics**: Navier-Stokes equations (linearized)
- **Heat Transfer**: Finite difference/element methods
- **Wave Propagation**:  Helmholtz equation

#### 8. **Data Science** 📊
- **Least Squares Fitting**: Overdetermined systems
- **Polynomial Interpolation**: Vandermonde systems
- **Statistical Inference**: Normal equations
- **Network Analysis**: Graph Laplacians

### Real-World Example: Power Grid

**Problem**: Calculate power flow in electrical grid

**Setup**:
- 1000 buses (nodes) in the network
- Need to solve for voltage at each bus
- Topology rarely changes, but loads change hourly

**Without LU**: 
- Solve 1000×1000 system 24 times per day
- Cost: 24 × (n³/3) ≈ 8 billion operations/day

**With LU**: 
- Decompose once per topology change (monthly)
- Solve 24 times per day:  24 × (2n²) ≈ 48 million operations/day
- **Speedup:  ~166×** for daily operations! 

---

## 📁 Implementation Structure

Each method folder contains:

```
Method Name/
├── README.md                 # Detailed theory and examples
├── method-name.cpp           # C++ implementation
├── input. txt                # Sample test cases
└── output.txt                # Expected outputs
```

### Common Features Across All Implementations

**Input/Output**:
- ✅ File-based I/O for batch processing
- ✅ Multiple test cases per file
- ✅ Formatted output with configurable precision

**Visualization**:
- ✅ Step-by-step matrix transformations
- ✅ Intermediate results display (togglable)
- ✅ Final solution verification

**Error Handling**: 
- ✅ Solution type detection (unique/none/infinite)
- ✅ Singular matrix detection
- ✅ Numerical stability checks (epsilon comparisons)
- ✅ Input validation

**Code Quality**:
- ✅ Well-commented code
- ✅ Consistent formatting
- ✅ Readable variable names
- ✅ Modular structure

---

## 🚀 Getting Started

### Prerequisites
- **C++ Compiler**: g++, clang++, or MSVC
- **C++ Standard**: C++11 or later
- **Text Editor/IDE**: VS Code, CLion, or any preferred editor
- **Basic Knowledge**: Linear algebra fundamentals

### Compilation

Navigate to any method folder and compile:

```bash
# Standard compilation
g++ -o solver method-name.cpp -std=c++11

# With optimization (recommended)
g++ -o solver method-name.cpp -std=c++17 -O2

# With warnings (for development)
g++ -o solver method-name.cpp -std=c++17 -O2 -Wall -Wextra
```

### Execution

```bash
# Run the compiled program
./solver

# On Windows
solver.exe
```

The program: 
1. Reads input from `input.txt`
2. Processes all test cases
3. Writes results to `output.txt`
4. Displays completion message

### Input Format

```
n                    # Number of equations
a₁₁ a₁₂ ... a₁ₙ b₁   # First equation coefficients and constant
a₂₁ a₂₂ ... a₂ₙ b₂   # Second equation coefficients and constant
...
aₙ₁ aₙ₂ ... aₙₙ bₙ     # nth equation coefficients and constant

# Multiple test cases can follow
```

**Example**:
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3

3
1 2 3 6
2 4 6 12
3 6 9 18
```

**Represents**:
```
Test Case 1:
2x + y - z = 8
-3x - y + 2z = -11
-2x + y + 2z = -3

Test Case 2:
x + 2y + 3z = 6
2x + 4y + 6z = 12
3x + 6y + 9z = 18
```

### Output Format

Each test case produces: 
1. **Input system display** (equations in readable form)
2. **Intermediate steps** (matrix transformations)
3. **Solution type** (Unique/None/Infinite)
4. **Final solution** (if unique)
5. **Verification** (Ax = b check)

### Customization Options

**Toggle Intermediate Output**:
```cpp
bool printIntermediate = true;  // Set to false for final results only
```

**Adjust Precision**:
```cpp
fout << fixed << setprecision(4);  // Change 4 to desired decimal places
```

**Modify Zero Threshold**:
```cpp
const double EPSILON = 1e-12;  // Adjust for your numerical requirements
```

---

## 📚 References

- **Numerical Methods for Engineers** by Chapra & Canale
- [Wikipedia: System of Linear Equations](https://en.wikipedia.org/wiki/System_of_linear_equations)

---

## 👨‍💻 Author

**Abir Hasan Arko**  
[![GitHub](https://img.shields.io/badge/GitHub-AbirHasanArko-181717?style=flat&logo=github)](https://github.com/AbirHasanArko)

---

<div align="center">

**[⬆ Back to Top](#solution-of-linear-equations)**

Part of the [Numerical Computing Suite](https://github.com/AbirHasanArko/Numerical-Computing-Suite)

</div>