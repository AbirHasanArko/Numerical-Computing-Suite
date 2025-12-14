# Interpolation and Approximation

[![C++](https://img.shields.io/badge/Language-C++-00599C?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)
[![Numerical Methods](https://img.shields.io/badge/Topic-Numerical%20Methods-FF6B6B?style=for-the-badge)](https://github.com/AbirHasanArko/Numerical-Computing-Suite)

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Overview of Interpolation Methods](#-overview-of-interpolation-methods)
- [Mathematical Foundation](#-mathematical-foundation)
  - [What is Interpolation?](#what-is-interpolation)
  - [What is Approximation?](#what-is-approximation)
  - [Polynomial Interpolation](#polynomial-interpolation)
  - [Fundamental Concepts](#fundamental-concepts)
- [Interpolation Methods](#-interpolation-methods)
  - [1. Newton's Forward Interpolation](#1-newtons-forward-interpolation)
  - [2. Newton's Backward Interpolation](#2-newtons-backward-interpolation)
  - [3. Newton's Divided Difference Interpolation](#3-newtons-divided-difference-interpolation)
- [Method Comparison](#-method-comparison)
- [Applications](#-applications)
- [Implementation Structure](#-implementation-structure)
- [Getting Started](#-getting-started)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

This section of the **Numerical Computing Suite** focuses on **interpolation and approximation** techniques for estimating unknown values between known data points.  Interpolation is a fundamental problem in numerical analysis that arises when we have discrete data points and need to estimate values at intermediate positions.

**The Core Problem**: Given a set of data points `(x₀, y₀), (x₁, y₁), ..., (xₙ, yₙ)`, find a function f(x) such that:
- f(xᵢ) = yᵢ for all known points (interpolation constraint)
- f(x) provides reasonable estimates for unknown x values

Interpolation appears in countless applications: 
- **Scientific Computing**: Fitting experimental data, signal processing
- **Computer Graphics**:  Smooth curve generation, animation
- **Engineering**: Data reconstruction, sensor calibration
- **Statistics**: Missing data estimation, trend analysis
- **Finance**: Option pricing, yield curve construction

This collection provides three powerful Newton interpolation methods, each optimized for different data characteristics and use cases.

---

## 🔍 Overview of Interpolation Methods

This repository implements three variants of **Newton's interpolation formula**:

| Method | Data Spacing | Best For | Complexity | Key Feature |
|--------|-------------|----------|------------|-------------|
| **Newton's Forward** | Equal spacing | Near beginning | O(n²) | Forward difference table |
| **Newton's Backward** | Equal spacing | Near end | O(n²) | Backward difference table |
| **Newton's Divided Difference** | Any spacing | General case | O(n²) | Works with non-uniform data |

---

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

## 🛠️ Interpolation Methods

## 1. Newton's Forward Interpolation

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-blue?style=for-the-badge)](./Newton's%20Forward%20Interpolation/)

### Theory

**Newton's Forward Interpolation** (also called **Newton's Forward Difference Formula**) is a polynomial interpolation method specifically designed for **equally spaced data points**. It's most accurate when interpolating near the **beginning** of the dataset.

### Historical Context

Named after Sir Isaac Newton (1642-1727), this method appeared in his work on finite differences in the 1670s. Newton recognized that equally spaced data allows for simpler formulas using difference operators, making hand calculations more practical before the computer age.

### Mathematical Foundation

**Assumptions**:
1. Data points are **equally spaced**:  xᵢ₊₁ - xᵢ = h (constant step size)
2. Interpolation point x is **near the beginning** of the data range
3. Function is reasonably smooth

**The Formula**: 

Given data points at x₀, x₁, x₂, ..., xₙ with equal spacing h, the interpolating polynomial is:

```
Pₙ(x) = y₀ + uΔy₀ + [u(u-1)/2! ]Δ²y₀ + [u(u-1)(u-2)/3!]Δ³y₀ + ... 
```

Where:
- **u = (x - x₀) / h** (normalized position)
- **Δʲy₀** = jth forward difference at x₀
- The formula uses differences at the **first point** (x₀)

**General Term**:
```
Term_k = [u(u-1)(u-2)...(u-k+1) / k!] × Δᵏy₀
```

#### Why "Forward" Difference?

The method is called "forward" because: 
- It uses the forward difference operator Δ
- Differences are computed moving forward from x₀
- Most accurate for interpolation near the start (where u is small)

**Forward Difference Table Structure**:

```
x     y       Δy      Δ²y     Δ³y     Δ⁴y
─────────────────────────────────────────
x₀    y₀      
             Δy₀
x₁    y₁             Δ²y₀
             Δy₁             Δ³y₀
x₂    y₂             Δ²y₁            Δ⁴y₀
             Δy₂             Δ³y₁
x₃    y₃             Δ²y₂
             Δy₃
x₄    y₄
```

**Key Observation**: All differences needed are in the **top row** or **diagonal from top**. 

### Detailed Algorithm

##### **Step 1: Validate Equal Spacing**

```
h = x₁ - x₀
For i = 1 to n-1:
    If |xᵢ₊₁ - xᵢ - h| > ε: 
        Error:  "Data points must be equally spaced"
```

##### **Step 2: Build Forward Difference Table**

```
Initialize table[n][n]
table[i][0] = yᵢ for all i  (first column = y values)

For j = 1 to n-1:           (difference order)
    For i = 0 to n-j-1:     (starting position)
        table[i][j] = table[i+1][j-1] - table[i][j-1]
```

**Example**: Data points (0, 1), (1, 2), (2, 5), (3, 10)

```
Initial: 
i   x   y=table[i][0]
0   0   1
1   1   2
2   2   5
3   3   10

After j=1 (first differences):
table[0][1] = table[1][0] - table[0][0] = 2 - 1 = 1
table[1][1] = table[2][0] - table[1][0] = 5 - 2 = 3
table[2][1] = table[3][0] - table[2][0] = 10 - 5 = 5

After j=2 (second differences):
table[0][2] = table[1][1] - table[0][1] = 3 - 1 = 2
table[1][2] = table[2][1] - table[1][1] = 5 - 3 = 2

After j=3 (third differences):
table[0][3] = table[1][2] - table[0][2] = 2 - 2 = 0

Final Table:
x   y    Δy   Δ²y  Δ³y
0   1    1    2    0
1   2    3    2
2   5    5
3   10
```

##### **Step 3: Interpolation**

For a given interpolation point x:

```
u = (x - x₀) / h

result = table[0][0]  (= y₀)
u_term = 1. 0
factorial = 1.0

For j = 1 to n-1:
    u_term = u_term × (u - (j-1))
    factorial = factorial × j
    term = (u_term / factorial) × table[0][j]
    result = result + term
    
    If |term| < ε:  // Early termination if contribution negligible
        break

Return result
```

#### Detailed Example Walkthrough

**Given**:  Approximate f(1.5) using the data: 

```
x:  1.0   1.5   2.0   2.5   3.0
y:  0.0   0.41  0.69  0.92  1.10
```

**Step 1**: Verify equal spacing
```
h = 1.5 - 1.0 = 0.5 ✓
All intervals = 0.5 ✓
```

**Step 2**: Build forward difference table

```
x    y      Δy      Δ²y      Δ³y      Δ⁴y
1.0  0.00   0.41    -0.13    0.02     -0.01
1.5  0.41   0.28    -0.11    0.01
2.0  0.69   0.23    -0.10
2.5  0.92   0.18
3.0  1.10
```

Calculations:
- Δy₀ = 0.41 - 0.00 = 0.41
- Δ²y₀ = 0.28 - 0.41 = -0.13
- Δ³y₀ = -0.11 - (-0.13) = 0.02
- Δ⁴y₀ = 0.01 - 0.02 = -0.01

**Step 3**: Interpolate at x = 1.5
```
u = (1.5 - 1.0) / 0.5 = 1.0

P(1.5) = y₀ + uΔy₀ + [u(u-1)/2!]Δ²y₀ + [u(u-1)(u-2)/3!]Δ³y₀ + ... 

Term 0: y₀ = 0.00
Term 1: 1.0 × 0.41 = 0.41
Term 2: [1.0(1.0-1)/2] × (-0.13) = 0.00
Term 3: [1.0(1.0-1)(1.0-2)/6] × 0.02 = 0.00
Term 4: [1.0(1.0-1)(1.0-2)(1.0-3)/24] × (-0.01) = 0.00

Result: P(1.5) = 0.00 + 0.41 = 0.41 ✓
```

(This matches the table value because 1.5 is an actual data point!)

**Step 4**: Interpolate at x = 1.25 (between data points)
```
u = (1.25 - 1.0) / 0.5 = 0.5

Term 0: 0.00
Term 1: 0.5 × 0.41 = 0.205
Term 2: [0.5(0.5-1)/2] × (-0.13) = [0.5×(-0.5)/2] × (-0.13) = 0.01625
Term 3: [0.5×(-0.5)×(-1. 5)/6] × 0.02 = 0.00125
Term 4: [0.5×(-0.5)×(-1.5)×(-2.5)/24] × (-0.01) ≈ -0.00039

Result: P(1.25) ≈ 0.00 + 0.205 + 0.016 + 0.001 ≈ 0.222
```

### When to Use Forward Interpolation

**Optimal Scenarios**:  
✅ Equally spaced data points  
✅ Interpolation near the **beginning** of the dataset (u ≈ 0 to 1)  
✅ Tabulated data (like mathematical tables)  
✅ Sequential data processing  

**Performance Characteristics**:
- **Most accurate**:  When 0 ≤ u ≤ 1 (first interval)
- **Good accuracy**: When 0 ≤ u ≤ 2 (first two intervals)
- **Degrading accuracy**: When u > 2 (far from start)

**Avoid When**:  
❌ Data points are not equally spaced (use divided differences)  
❌ Interpolating near the end of dataset (use backward interpolation)  
❌ High-degree polynomial needed (risk of Runge's phenomenon)  

### Complexity Analysis

**Time Complexity**: 

1. **Building Difference Table** (one-time):
   - Outer loop: j = 1 to n-1
   - Inner loop: i = 0 to n-j-1
   - Total operations:  Σⱼ₌₁ⁿ⁻¹(n-j) = n(n-1)/2
   - **Complexity**: O(n²)

2. **Single Interpolation** (using pre-built table):
   - Loop: k = 1 to n-1
   - Each iteration: O(1) operations
   - **Complexity**:  O(n)

3. **Total for k Interpolations**:
   - Build once: O(n²)
   - k interpolations: O(kn)
   - **Overall**: O(n² + kn)

**Space Complexity**:
- Difference table: n×n matrix (though only upper triangle used)
- **Space**: O(n²)
- Can be optimized to O(n) if only storing diagonal

**Operation Counts** (for single interpolation with n points):
- **Multiplications**: ~n² (table) + ~n² (interpolation) ≈ 2n²
- **Additions/Subtractions**: ~n²/2 (table) + ~n (interpolation)
- **Divisions**: ~n (computing u terms and factorials)

### Numerical Considerations

**Sources of Error**: 

1. **Truncation Error**:
   - From stopping at finite degree
   - Magnitude:  O(hⁿ⁺¹) where h is step size
   - Solution: Use more data points or smaller h

2. **Round-off Error**:
   - From finite precision arithmetic
   - Accumulates in difference table
   - Solution: Use double precision, check for constant high-order differences

3. **Runge's Phenomenon**:
   - High-degree polynomials oscillate
   - Especially problematic with many points
   - Solution:  Limit polynomial degree or use piecewise interpolation

**Stability Tips**:  
✅ Verify equal spacing with tolerance check  
✅ Check for constant high-order differences (indicates polynomial degree)  
✅ Use double precision for all calculations  
✅ Validate interpolation point is within reasonable range  
✅ Limit polynomial degree to avoid oscillations  

#### Advantages & Limitations

**Advantages**:  
✅ **Simple and intuitive** formula  
✅ **Efficient** for equally spaced data  
✅ **Reusable table** for multiple interpolations  
✅ **Pattern recognition** - constant differences reveal polynomial degree  
✅ **Natural fit** for tabulated data  
✅ **Easy error estimation** from high-order terms  
✅ **Historical tables** already in this format  

**Limitations**:  
❌ **Requires equal spacing** - strict constraint  
❌ **Best near start** - accuracy degrades away from x₀  
❌ **High-degree issues** - oscillations with many points  
❌ **Memory intensive** - O(n²) storage  
❌ **Not general-purpose** - specific use case  

---

## 2. Newton's Backward Interpolation

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-green?style=for-the-badge)](./Newton's%20Backward%20Interpolation/)

### Theory

**Newton's Backward Interpolation** (also called **Newton's Backward Difference Formula**) is the mirror image of forward interpolation. It's specifically designed for **equally spaced data points** and is most accurate when interpolating near the **end** of the dataset.

### Historical Context

Developed alongside the forward formula, backward interpolation became essential for astronomical and navigational tables where extrapolation beyond the last known value was common.  By working backward from the end, it provides better accuracy for such calculations.

### Mathematical Foundation

**Assumptions**:
1. Data points are **equally spaced**: xᵢ₊₁ - xᵢ = h (constant step size)
2. Interpolation point x is **near the end** of the data range
3. Function is reasonably smooth

**The Formula**:

Given data points at x₀, x₁, x₂, ..., xₙ with equal spacing h, the interpolating polynomial is: 

```
Pₙ(x) = yₙ + u∇yₙ + [u(u+1)/2!]∇²yₙ + [u(u+1)(u+2)/3!]∇³yₙ + ...
```

Where:
- **u = (x - xₙ) / h** (normalized position from the end)
- **∇ʲyₙ** = jth backward difference at xₙ (last point)
- The formula uses differences at the **last point** (xₙ)

**General Term**:
```
Term_k = [u(u+1)(u+2)...(u+k-1) / k!] × ∇ᵏyₙ
```

**Key Difference from Forward**:
- Forward uses:  u(u-1)(u-2)... and Δʲy₀ (from start)
- Backward uses: u(u+1)(u+2)... and ∇ʲyₙ (from end)

#### Why "Backward" Difference? 

The method is called "backward" because:
- It uses the backward difference operator ∇
- Differences are computed moving backward from xₙ
- Most accurate for interpolation near the end (where u is small in magnitude)

**Backward Difference Table Structure**:

```
x     y       ∇y      ∇²y     ∇³y     ∇⁴y
─────────────────────────────────────────
x₀    y₀      
              ∇y₁
x₁    y₁             ∇²y₂
              ∇y₂             ∇³y₃
x₂    y₂             ∇²y₃             ∇⁴y₄
              ∇y₃             ∇³y₄
x₃    y₃             ∇²y₄
              ∇y₄
x₄    y₄    ←───────────────── Bottom row
```

**Key Observation**: All differences needed are in the **bottom row** or **diagonal to bottom**.

### Detailed Algorithm

##### **Step 1: Validate Equal Spacing**

```
h = x₁ - x₀
For i = 1 to n-1:
    If |xᵢ₊₁ - xᵢ - h| > ε: 
        Error: "Data points must be equally spaced"
```

##### **Step 2: Build Backward Difference Table**

```
Initialize table[n][n]
table[i][0] = yᵢ for all i  (first column = y values)

For j = 1 to n-1:           (difference order)
    For i = n-1 down to j:  (starting from bottom)
        table[i][j] = table[i][j-1] - table[i-1][j-1]
```

**Key Difference from Forward**:  The inner loop goes **bottom-to-top** (i = n-1 down to j).

**Example**: Data points (0, 1), (1, 2), (2, 5), (3, 10)

```
Initial:
i   x   y=table[i][0]
0   0   1
1   1   2
2   2   5
3   3   10

After j=1 (first backward differences):
table[3][1] = table[3][0] - table[2][0] = 10 - 5 = 5
table[2][1] = table[2][0] - table[1][0] = 5 - 2 = 3
table[1][1] = table[1][0] - table[0][0] = 2 - 1 = 1

After j=2 (second backward differences):
table[3][2] = table[3][1] - table[2][1] = 5 - 3 = 2
table[2][2] = table[2][1] - table[1][1] = 3 - 1 = 2

After j=3 (third backward differences):
table[3][3] = table[3][2] - table[2][2] = 2 - 2 = 0

Final Backward Difference Table:
x   y    ∇y   ∇²y  ∇³y
0   1    
1   2    1
2   5    3    2
3   10   5    2    0  ← bottom row contains all needed differences
```

##### **Step 3: Interpolation**

For a given interpolation point x near the end:

```
u = (x - xₙ) / h  // Note: typically negative for interpolation

result = table[n-1][0]  (= yₙ, last y value)
u_term = 1.0
factorial = 1.0

For j = 1 to n-1:
    u_term = u_term × (u + (j-1))  // Note: u+0, u+1, u+2, ...
    factorial = factorial × j
    term = (u_term / factorial) × table[n-1][j]  // Use bottom row
    result = result + term
    
    If |term| < ε: 
        break

Return result
```

**Critical Detail**: u_term multiplies by (u + j-1), not (u - j+1) like forward! 

#### Detailed Example Walkthrough

**Given**: Approximate f(2.75) using the data:

```
x:  1.0   1.5   2.0   2.5   3.0
y:  0.0   0.41  0.69  0.92  1.10
```

**Step 1**: Verify equal spacing
```
h = 0.5 ✓
```

**Step 2**: Build backward difference table

```
x    y      ∇y      ∇²y      ∇³y      ∇⁴y
1.0  0.00   
1.5  0.41   0.41   
2.0  0.69   0.28    -0.13
2.5  0.92   0.23    -0.05    0.08
3.0  1.10   0.18    -0.05    0.00     -0.08  ← bottom row
```

Backward differences (calculated bottom-up):
- ∇y₄ = y₄ - y₃ = 1.10 - 0.92 = 0.18
- ∇²y₄ = ∇y₄ - ∇y₃ = 0.18 - 0.23 = -0.05
- ∇³y₄ = ∇²y₄ - ∇²y₃ = -0.05 - (-0.05) = 0.00
- ∇⁴y₄ = ∇³y₄ - ∇³y₃ = 0.00 - 0.08 = -0.08

**Step 3**: Interpolate at x = 2.75
```
n = 4 (index of last point is 4)
xₙ = 3.0
u = (2.75 - 3.0) / 0.5 = -0.5

P(2.75) = yₙ + u∇yₙ + [u(u+1)/2!]∇²yₙ + [u(u+1)(u+2)/3!]∇³yₙ + ...

Term 0: y₄ = 1.10
Term 1: (-0.5) × 0.18 = -0.09
Term 2: [(-0.5)(0.5)/2] × (-0.05) = [-0.25/2] × (-0.05) = 0.00625
Term 3: [(-0.5)(0.5)(1.5)/6] × 0.00 = 0.00
Term 4: [(-0.5)(0.5)(1.5)(2.5)/24] × (-0.08) ≈ 0.0039

Result: P(2.75) ≈ 1.10 - 0.09 + 0.006 + 0.004 ≈ 1.020
```

**Verification**: This is between y₃ = 0.92 and y₄ = 1.10, which makes sense!  ✓

### Forward vs Backward:  Direct Comparison

**Same Polynomial, Different Representation**: 

Both methods produce the **same interpolating polynomial**, just expressed differently! 

| Aspect | Forward | Backward |
|--------|---------|----------|
| **Reference Point** | First point (x₀) | Last point (xₙ) |
| **Parameter u** | (x - x₀) / h | (x - xₙ) / h |
| **Typical u Range** | 0 to n | -n to 0 |
| **Difference Operator** | Δ (forward) | ∇ (backward) |
| **Table Row Used** | Top row | Bottom row |
| **Terms** | u(u-1)(u-2)... | u(u+1)(u+2)... |
| **Best Accuracy** | Near start (small u) | Near end (u ≈ 0) |
| **Use Case** | Interpolation at beginning | Interpolation at end |

**Example**:  For the same data and interpolation point in the middle, both give the same result (within numerical precision).

### When to Use Backward Interpolation

**Optimal Scenarios**:  
✅ Equally spaced data points  
✅ Interpolation near the **end** of the dataset  
✅ **Extrapolation** slightly beyond last point (with caution)  
✅ Time series where most recent data is most relevant  
✅ Sequential processing from end to start  

**Performance Characteristics**:
- **Most accurate**: When -1 ≤ u ≤ 0 (last interval)
- **Good accuracy**:  When -2 ≤ u ≤ 0 (last two intervals)
- **Degrading accuracy**: When u < -2 (far from end)

**Real-World Example**:
Stock price prediction:  If you have daily prices and want to estimate tomorrow's price (extrapolation beyond last point), backward interpolation is more appropriate than forward.

### Complexity Analysis

**Time Complexity**: 
- **Building Table**: O(n²) - same as forward
- **Single Interpolation**: O(n) - same as forward
- **k Interpolations**: O(n² + kn)

**Space Complexity**:
- O(n²) for full table
- Can optimize to O(n) storing only bottom diagonal

**Operation Counts**:  Identical to forward interpolation

#### Numerical Considerations

**Error Sources**:  Same as forward interpolation
- Truncation error from finite degree
- Round-off error in arithmetic
- Runge's phenomenon for high degrees

**Extrapolation Warning**:
```
⚠️  Extrapolating beyond xₙ (u > 0) can be very inaccurate! 
    Use with extreme caution and check reasonableness. 
```

**Why Extrapolation is Risky**:
- No data to constrain polynomial beyond last point
- High-degree terms dominate
- Small errors in data cause large errors in extrapolation

#### Advantages & Limitations

**Advantages**:  
✅ **Mirror of forward** - same theoretical properties  
✅ **Better for end points** - accuracy where forward is weak  
✅ **Complements forward** - choose based on interpolation location  
✅ **Same efficiency** - O(n²) table, O(n) interpolation  
✅ **Historical importance** - essential for navigation tables  

**Limitations**:  
❌ **Same restrictions as forward** - needs equal spacing  
❌ **Not general-purpose** - specific to end-region interpolation  
❌ **Extrapolation dangers** - easily misused  
❌ **Memory intensive** - O(n²) storage  

---

## 3. Newton's Divided Difference Interpolation

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-orange?style=for-the-badge)](./Newton's%20Divided%20Difference%20Interpolation/)

### Theory

**Newton's Divided Difference Interpolation** is the most **general form** of Newton's interpolation method. Unlike forward and backward interpolation which require equally spaced data, this method works with **arbitrarily spaced data points**, making it the most versatile and widely applicable interpolation technique.

### Historical Context

This is Newton's original formulation of polynomial interpolation from his work in the 1670s-1680s. The method predates the specialized forward/backward formulas, which were later simplified for the common case of equal spacing.  The divided difference approach represents the foundational theory underlying all Newton interpolation methods.

### Mathematical Foundation

**Key Advantage**: **No assumption about data spacing** - works with any distinct x-values! 

**Assumptions**:
1. Data points have **distinct x-values** (no duplicates)
2. Points can be in **any order** (though sorted is conventional)
3. Spacing can be **completely arbitrary**
4. Function is reasonably smooth

**The Formula**:

Given data points (x₀, y₀), (x₁, y₁), ..., (xₙ, yₙ), the interpolating polynomial is:

```
Pₙ(x) = f[x₀] + f[x₀,x₁](x-x₀) + f[x₀,x₁,x₂](x-x₀)(x-x₁) + ... 
        + f[x₀,... ,xₙ](x-x₀)(x-x₁).. .(x-xₙ₋₁)
```

Where **f[x₀, x₁, ..., xₖ]** is the **kth divided difference**. 

**Compact Form**:
```
Pₙ(x) = Σₖ₌₀ⁿ f[x₀,... ,xₖ] ∏ⱼ₌₀ᵏ⁻¹(x - xⱼ)
```

### Divided Differences:  The Core Concept

**Divided differences** are a generalization of derivatives and differences that work for non-uniform spacing. 

**Recursive Definition**: 

**Order 0** (zeroth divided difference):
```
f[xᵢ] = f(xᵢ) = yᵢ
```
Just the function value! 

**Order 1** (first divided difference):
```
f[xᵢ, xᵢ₊₁] = (f[xᵢ₊₁] - f[xᵢ]) / (xᵢ₊₁ - xᵢ)
            = (yᵢ₊₁ - yᵢ) / (xᵢ₊₁ - xᵢ)
```
This is the **slope** between two points!

**Order 2** (second divided difference):
```
f[xᵢ, xᵢ₊₁, xᵢ₊₂] = (f[xᵢ₊₁,xᵢ₊₂] - f[xᵢ,xᵢ₊₁]) / (xᵢ₊₂ - xᵢ)
```

**General Order k**:
```
f[xᵢ, xᵢ₊₁, .. ., xᵢ₊ₖ] = (f[xᵢ₊₁,... ,xᵢ₊ₖ] - f[xᵢ,... ,xᵢ₊ₖ₋₁]) / (xᵢ₊ₖ - xᵢ)
```

**Key Properties**:
1. **Symmetric**: f[x₀, x₁, x₂] = f[x₁, x₀, x₂] = f[x₂, x₁, x₀] (order doesn't matter!)
2. **Connection to derivatives**: If spacing → 0, f[x,x+h] → f'(x)
3. **Polynomial detection**: For polynomial of degree n, (n+1)th divided differences are constant
4. **Reduces to forward differences**:  When h is constant, f[x₀, x₁] = Δy₀/h

#### Why "Divided" Difference?

Called "divided" because we **divide by the x-interval** (xᵢ₊ₖ - xᵢ), making it a **rate of change** that generalizes to non-uniform spacing.

**Intuition**:
- First divided difference: average rate of change (slope)
- Second divided difference:  rate of change of rate of change (curvature)
- Higher orders:  increasingly subtle shape information

#### Divided Difference Table Structure

Unlike forward/backward tables, the divided difference table has a different structure:

```
x       f[xᵢ]    f[xᵢ,xᵢ₊₁]   f[xᵢ,xᵢ₊₁,xᵢ₊₂]   f[x₀... x₃]   f[x₀...x₄]
──────────────────────────────────────────────────────────────────────────
x₀      y₀       f[x₀,x₁]      f[x₀,x₁,x₂]       f[x₀... x₃]   f[x₀...x₄]
x₁      y₁       f[x₁,x₂]      f[x₁,x₂,x₃]       f[x₁...x₄]
x₂      y₂       f[x₂,x₃]      f[x₂,x₃,x₄]
x₃      y₃       f[x₃,x₄]
x₄      y₄
```

**Key Observation**: The **top row** contains all coefficients needed for the interpolating polynomial! 

### Detailed Algorithm

##### **Step 1: Validate Input**

```
For i = 0 to n-1:
    For j = i+1 to n:
        If |xᵢ - xⱼ| < ε:
            Error:  "Duplicate x-values detected"
```

**Why This Matters**:  Duplicate x-values cause **division by zero** in divided difference formula!

##### **Step 2: Build Divided Difference Table**

```
Initialize table[n][n]

// Column 0: function values
For i = 0 to n-1:
    table[i][0] = yᵢ

// Build subsequent columns
For j = 1 to n-1:              // Column (difference order)
    For i = 0 to n-j-1:        // Row (starting position)
        numerator = table[i+1][j-1] - table[i][j-1]
        denominator = x[i+j] - x[i]
        table[i][j] = numerator / denominator
```

**Example**: Data points (1, 1), (2, 8), (4, 64), (5, 125)

```
Step-by-step construction: 

Column 0 (f[xᵢ]):
table[0][0] = 1
table[1][0] = 8
table[2][0] = 64
table[3][0] = 125

Column 1 (f[xᵢ, xᵢ₊₁]):
table[0][1] = (8 - 1) / (2 - 1) = 7
table[1][1] = (64 - 8) / (4 - 2) = 28
table[2][1] = (125 - 64) / (5 - 4) = 61

Column 2 (f[xᵢ, xᵢ₊₁, xᵢ₊₂]):
table[0][2] = (28 - 7) / (4 - 1) = 7
table[1][2] = (61 - 28) / (5 - 2) = 11

Column 3 (f[x₀, x₁, x₂, x₃]):
table[0][3] = (11 - 7) / (5 - 1) = 1

Final Table:
x    f[xᵢ]   f[.,.]   f[. ,. ,.] f[.,.,.,. ] 
1    1       7        7        1
2    8       28       11
4    64      61
5    125
```

**Pattern Recognition**:  These are coefficients of x³ (since constant 4th divided difference doesn't exist here with 4 points).

##### **Step 3: Interpolation Using Nested Multiplication (Horner's Method)**

For a given interpolation point x, evaluate Pₙ(x) efficiently:

```
// Horner's method (most efficient)
result = table[n-1][n-1]  // Start with highest-order term
For i = n-2 down to 0:
    result = result × (x - x[i]) + table[i][i]

Return result
```

**Alternative (Direct Evaluation)**:
```
result = table[0][0]
term = 1.0

For k = 1 to n-1:
    term = term × (x - x[k-1])
    result = result + table[0][k] × term

Return result
```

Both methods give the same result; Horner's is slightly more efficient. 

#### Detailed Example Walkthrough

**Given**: Approximate f(3) using data:  (0, 1), (0.7, 2.014), (1.3, 3.669), (2.0, 7.389)

These represent f(x) = e^x evaluated at unequally spaced points.

**Step 1**:  Validate - no duplicate x-values ✓

**Step 2**:  Build divided difference table

```
x     f[xᵢ]    f[xᵢ,xᵢ₊₁]   f[xᵢ,xᵢ₊₁,xᵢ₊₂]   f[x₀...x₃]
────────────────────────────────────────────────────────────
0. 0   1.000    1.449        1.008              0.479
0.7   2.014    2.758        1.966
1.3   3.669    5.314
2.0   7.389
```

Calculations:
```
f[x₀,x₁] = (2.014 - 1.000) / (0.7 - 0.0) = 1.014 / 0.7 = 1.449
f[x₁,x₂] = (3.669 - 2.014) / (1.3 - 0.7) = 1.655 / 0.6 = 2.758
f[x₂,x₃] = (7.389 - 3.669) / (2.0 - 1.3) = 3.720 / 0.7 = 5.314

f[x₀,x₁,x₂] = (2.758 - 1.449) / (1.3 - 0.0) = 1.309 / 1.3 = 1.008
f[x₁,x₂,x₃] = (5.314 - 2.758) / (2.0 - 0.7) = 2.556 / 1.3 = 1.966

f[x₀,x₁,x₂,x₃] = (1.966 - 1.008) / (2.0 - 0.0) = 0.958 / 2.0 = 0.479
```

**Step 3**: Interpolate at x = 3

Using direct evaluation:
```
P₃(3) = f[x₀] 
      + f[x₀,x₁](3-x₀) 
      + f[x₀,x₁,x₂](3-x₀)(3-x₁) 
      + f[x₀,x₁,x₂,x₃](3-x₀)(3-x₁)(3-x₂)

     = 1.000 
     + 1.449 × (3 - 0) 
     + 1.008 × (3 - 0) × (3 - 0.7) 
     + 0.479 × (3 - 0) × (3 - 0.7) × (3 - 1.3)

     = 1.000 
     + 1.449 × 3 
     + 1.008 × 3 × 2.3 
     + 0.479 × 3 × 2.3 × 1.7

     = 1.000 + 4.347 + 6.955 + 5.606
     = 17.908
```

**Verification**: The actual value is e³ ≈ 20.086. Our estimate 17.908 is reasonable given we're extrapolating! 

#### Connection to Other Forms

**Relation to Lagrange Interpolation**:
- Newton and Lagrange produce the **same polynomial**
- Newton:  Builds incrementally, easy to add points
- Lagrange: Direct formula, harder to add points

**Relation to Forward/Backward Differences**:
When spacing is equal (h = constant):
```
f[xᵢ, xᵢ₊₁] = Δyᵢ / h
f[xᵢ, xᵢ₊₁, xᵢ₊₂] = Δ²yᵢ / (2! h²)
f[x₀,... ,xₖ] = Δᵏy₀ / (k!hᵏ)
```

So forward/backward interpolation are **special cases** of divided differences! 

#### Key Advantages Over Forward/Backward

**1. Flexibility in Data**:
```
Forward/Backward:   Must have x₀, x₀+h, x₀+2h, ...  (rigid)
Divided Difference: Any x₀, x₁, x₂, ... (flexible) ✓
```

**2. Incremental Updates**:
Adding a new point (xₙ₊₁, yₙ₊₁):
- Forward/Backward: Must rebuild entire table if spacing changes
- Divided Difference:  Just add one more column ✓

**3. Real-World Data**:
Most real data doesn't have perfect equal spacing!

### When to Use Divided Differences

**Optimal Scenarios**:  
✅ **Non-uniformly spaced data** (most common case)  
✅ **Experimental data** with irregular sampling  
✅ **General-purpose interpolation** - works everywhere  
✅ **Adding points incrementally** - table extends easily  
✅ **Adaptive algorithms** - can choose optimal points  
✅ **Legacy data** with missing or irregular entries  

**Real-World Examples**:
- Temperature readings at irregular times
- Stock prices (markets closed on weekends/holidays)
- Scientific measurements at varying intervals
- Historical data with gaps

**Still Use Forward/Backward When**:  
❌ Data is perfectly equally spaced (forward/backward slightly simpler)  
❌ Pedagogical purposes (easier to understand differences)  
❌ Historical compatibility (old tables in those formats)  

### Complexity Analysis

**Time Complexity**: 

1. **Building Table**:
   ```
   For j = 1 to n-1:      // n-1 iterations
       For i = 0 to n-j-1: // decreasing iterations
           // O(1) operation
   Total:  Σⱼ₌₁ⁿ⁻¹(n-j) = n(n-1)/2
   ```
   **Complexity**:  O(n²)

2. **Single Interpolation**:
   - Direct method: O(n) operations
   - Horner's method: O(n) operations (slightly fewer multiplications)

3. **k Interpolations**:  O(n² + kn)

**Space Complexity**:
- Full table: O(n²)
- Optimized (only top row needed): O(n)

**Operation Counts** (single interpolation, n points):
- **Divisions**: ~n²/2 (building table)
- **Multiplications**: ~n²/2 (table) + ~n (evaluation)
- **Additions/Subtractions**: ~n²/2 (table) + ~n (evaluation)

### Numerical Considerations

**Error Sources**: 

1. **Subtraction in Numerator**:
   - When f[xᵢ₊₁,ⱼ₋₁] ≈ f[xᵢ,ⱼ₋₁], subtractive cancellation occurs
   - Loss of significant digits
   - **Mitigation**: Use higher precision

2. **Small Denominators**:
   - When xᵢ₊ⱼ - xᵢ is small (nearly duplicate points)
   - Amplifies numerator errors
   - **Mitigation**:  Remove near-duplicate points, validate spacing

3. **Condition Number**:
   - Poorly spaced points increase condition number
   - **Best**:  Chebyshev points (cosine distribution)
   - **Worst**: Equally spaced endpoints (Runge's phenomenon)

**Stability Tips**:  
✅ Check for duplicate x-values before building table  
✅ Use double precision (at least)  
✅ Avoid extreme non-uniformity in spacing  
✅ Monitor high-order divided differences (should decrease)  
✅ Consider Chebyshev points for better conditioning  

### Error Analysis

**Interpolation Error Formula**:
```
E(x) = f(x) - Pₙ(x) = f[x₀,x₁,...,xₙ,x] × ∏ᵢ₌₀ⁿ(x - xᵢ)
```

Where f[x₀,... ,xₙ,x] is the (n+1)th divided difference.

**Practical Error Bound**:
If |f⁽ⁿ⁺¹⁾(ξ)| ≤ M for some ξ in [min xᵢ, max xᵢ]: 
```
|E(x)| ≤ (M / (n+1)!) × |∏ᵢ₌₀ⁿ(x - xᵢ)|
```

**Implications**:
- Error depends on smoothness of f (via f⁽ⁿ⁺¹⁾)
- Product term is zero at data points (exact interpolation)
- Error smallest near cluster of data points
- Error largest in gaps between data points

### Advantages & Limitations

**Advantages**:  
✅ **Most general** Newton method - no spacing restrictions  
✅ **Real-world ready** - works with actual data  
✅ **Incremental updates** - easy to add new points  
✅ **Flexible** - handles any distinct x-values  
✅ **Same polynomial** as Lagrange but easier to compute  
✅ **Error estimation** - from high-order terms  
✅ **Industry standard** - most numerical libraries use this  
✅ **Numerically stable** with proper implementation  

**Limitations**:  
❌ **More complex** than forward/backward (but not much)  
❌ **Duplicate detection needed** - must validate input  
❌ **Slightly more operations** than specialized equal-spacing methods  
❌ **Runge's phenomenon** still possible with many points  
❌ **Not optimal for derivatives** - finite differences better for that  

**Bottom Line**:  This is the **go-to method** for general interpolation! 

---

## 📊 Method Comparison

### Comprehensive Comparison Table

| Aspect | Forward | Backward | Divided Difference |
|--------|---------|----------|-------------------|
| **Spacing Requirement** | Equal (strict) | Equal (strict) | Any (flexible) ✓ |
| **Reference Point** | First (x₀) | Last (xₙ) | All points |
| **Best Accuracy Region** | Near start | Near end | Everywhere |
| **Difference Operator** | Δ (forward) | ∇ (backward) | Divided (/) |
| **Table Structure** | Top row | Bottom row | Top row |
| **Table Complexity** | O(n²) | O(n²) | O(n²) |
| **Interpolation Complexity** | O(n) | O(n) | O(n) |
| **Memory** | O(n²) | O(n²) | O(n²) |
| **Incremental Updates** | Hard | Hard | Easy ✓ |
| **General Purpose** | No | No | Yes ✓ |
| **Historical Usage** | Math tables | Navigation | Modern computing |
| **Implementation** | Simple | Simple | Moderate |

### When to Use Each Method? 

**Decision Tree**: 

```
Are data points equally spaced?
├─ YES → Where do you need to interpolate?
│        ├─ Near beginning → Newton's Forward
│        ├─ Near end → Newton's Backward
│        └─ Throughout → Any method works, divided difference most flexible
│
└─ NO → Newton's Divided Difference (only option)
```

### Accuracy Comparison

For the same data and interpolation point, all three methods produce the **same polynomial** (when applicable). The difference is in: 

1. **Numerical stability**:  Divided difference can be slightly less stable with poor spacing
2. **Computation efficiency**: Forward/backward slightly faster with equal spacing
3. **Applicability**: Only divided difference works universally

**Example**: Interpolating sin(x) at x = 0.25 using 5 points

```
Method                    Result          Error        Speed
Forward (0 to 1)         0.247404        3.0e-6       Fast
Backward (0 to 1)        0.247404        3.0e-6       Fast
Divided Diff (0 to 1)    0.247404        3.0e-6       Fast
Divided Diff (irregular) 0.247398        9.2e-6       Fast
```

All are accurate, but divided difference handles both cases! 

### Computational Cost (n = 10 points)

| Operation | Forward | Backward | Divided Diff |
|-----------|---------|----------|--------------|
| Table Build | 45 ops | 45 ops | 45 ops |
| Interpolate (1x) | 10 ops | 10 ops | 10 ops |
| Add New Point | Rebuild (45) | Rebuild (45) | Add column (10) ✓ |

---

## 🎯 Applications

### Scientific and Engineering Applications

#### 1. **Experimental Data Analysis** 🔬

**Problem**: Laboratory measurements at irregular time intervals

**Scenario**:
- Chemical reaction sampled at t = 0, 0.5, 1.2, 2.1, 3.5, 5.0 minutes
- Need concentration at t = 1.5 minutes

**Method**: Newton's Divided Difference (irregular spacing)

**Why Interpolation?**:  More reliable than fitting a complex model when underlying physics is partially known. 

#### 2. **Signal Processing** 📡

**Problem**: Reconstruct continuous signal from discrete samples

**Applications**:
- Audio upsampling (increase sample rate)
- Image scaling (enlarge photos)
- Video frame interpolation (slow motion)

**Method**: Typically divided difference or specialized methods (cubic splines)

**Example**: Convert 44.1 kHz audio to 48 kHz for video synchronization. 

#### 3. **Weather Prediction** 🌦️

**Problem**: Estimate temperature/pressure between weather stations

**Scenario**:
- Stations at irregular geographic positions
- Need values for grid-based numerical models

**Method**:  Divided difference for spatial interpolation

**Real Impact**: Improves accuracy of weather forecasting models! 

#### 4. **Engineering Tables** 📐

**Historical Application**: Before calculators and computers

**Example**: Steam tables
- Pressure, temperature, enthalpy at specific points
- Engineers interpolated for intermediate values

**Method**: Forward/backward (tables were equally spaced for this reason!)

**Modern Use**: Still used for quick estimates! 

### Computer Science Applications

#### 5. **Computer Graphics** 🎨

**Curve Generation**:
- User clicks 5 points
- Computer draws smooth curve through them
- Interpolation creates intermediate points

**Animation**:
- Keyframe at t = 0, 1, 3, 5 seconds
- Interpolate positions at t = 0.1, 0.2, ... 
- Creates smooth motion

**Font Rendering**:
- Outline defined by control points
- Interpolation generates smooth curves
- TrueType fonts use polynomial curves

#### 6. **Game Development** 🎮

**Camera Paths**:
- Define camera positions at key moments
- Interpolate smooth path between them
- Player experiences cinematic movement

**Particle Systems**:
- Particle properties (size, color, speed) over time
- Interpolate between keyframes
- Creates complex effects efficiently

#### 7. **Data Visualization** 📈

**Smooth Plotting**:
- Plot 10 data points
- Draw smooth curve through them
- Makes trends more visible

**Missing Data**:
- Sensor failed at some readings
- Interpolate to fill gaps
- Enables continuous analysis

### Financial Applications

#### 8. **Option Pricing** 💹

**Problem**: Option price depends on volatility surface

**Data**: Known prices at specific strikes/maturities
**Need**: Price for arbitrary strike/maturity
**Method**: 2D interpolation (combination of 1D interpolations)

**Impact**: Billions of dollars in daily trading rely on accurate interpolation!

#### 9. **Yield Curve Construction** 📊

**Problem**: Interest rates known at specific maturities

**Example**:
- 1-month:  2.5%
- 3-month: 2.7%
- 6-month:  2.9%
- 1-year:  3.2%

**Need**: Rate for 4-month loan

**Method**: Divided difference interpolation

**Use**:  Pricing bonds, derivatives, risk management

### Medical Applications

#### 10. **Medical Imaging** 🏥

**CT/MRI Scans**:
- Slices captured at specific intervals (e.g., every 5mm)
- Interpolate to estimate tissue density between slices
- Enables 3D reconstruction

**Radiation Treatment**:
- Plan radiation dose at specific points
- Interpolate to ensure proper dose distribution
- Critical for targeting tumors while sparing healthy tissue

**Prosthetic Design**:
- Measure body contours at specific locations
- Interpolate to create smooth 3D model
- Enables custom-fitted prosthetics

### Aerospace Applications

#### 11. **Trajectory Planning** 🚀

**Problem**: Spacecraft position at specific times

**Given**: Position at t = 0, 100, 200, 300 seconds
**Need**: Position at t = 150 seconds for mid-course correction

**Method**:  Divided difference (time intervals may vary)

**Impact**: Essential for navigation and rendezvous operations

#### 12. **Aerodynamic Tables** ✈️

**Problem**: Lift/drag coefficients measured at specific angles

**Data**: Wind tunnel tests at α = 0°, 5°, 10°, 15°, 20°
**Need**:  Coefficient at α = 12° for flight simulation

**Method**: Newton interpolation

**Use**: Aircraft design, flight simulators, autopilot systems

### Environmental Science

#### 13. **Climate Modeling** 🌍

**Temperature Reconstruction**:
- Historical data at irregular intervals
- Ice core samples, tree rings (non-uniform time spacing)
- Interpolate to create continuous temperature record

**Pollution Mapping**:
- Sensors at scattered locations
- Interpolate to estimate pollution levels city-wide
- Informs public health decisions

---

## 📁 Implementation Structure

Each method folder contains: 

```
Method Name/
├── README.md                              # Detailed theory and examples
├── newtons-[method]-interpolation.cpp    # C++ implementation
├── input1.txt                            # Basic test case
├── input2.txt                            # Error analysis case
├── output1.txt                           # Expected output 1
├── output2.txt                           # Expected output 2
```

### Common Features Across All Implementations

**Input/Output**:
- ✅ File-based I/O for batch processing
- ✅ Multiple test cases per file
- ✅ Formatted output with configurable precision
- ✅ Dual output (console and file)

**Visualization**:
- ✅ Data points table display
- ✅ Complete difference/divided difference table
- ✅ Step-by-step interpolation results
- ✅ Optional error analysis with additional points

**Error Handling**:
- ✅ Equal spacing validation (forward/backward)
- ✅ Duplicate x-value detection (divided difference)
- ✅ Extrapolation warnings
- ✅ Input validation and meaningful error messages

**Code Quality**:
- ✅ Well-commented code
- ✅ Consistent formatting
- ✅ Readable variable names
- ✅ Modular structure with utility functions

---

## 🚀 Getting Started

### Prerequisites
- **C++ Compiler**: g++, clang++, or MSVC
- **C++ Standard**: C++11 or later
- **Text Editor/IDE**: VS Code, CLion, or any preferred editor
- **Basic Knowledge**: Understanding of polynomials and basic calculus

### Compilation

Navigate to any method folder and compile:

```bash
# Standard compilation
g++ -o interpolate newtons-[method]-interpolation. cpp -std=c++11

# With optimization (recommended)
g++ -o interpolate newtons-[method]-interpolation. cpp -std=c++17 -O2

# With warnings (for development)
g++ -o interpolate newtons-[method]-interpolation.cpp -std=c++17 -O2 -Wall -Wextra
```

### Execution

```bash
# Run the compiled program
./interpolate

# On Windows
interpolate.exe
```

The program: 
1. Reads input from `input1.txt` or `input2.txt`
2. Processes all interpolation requests
3. Writes results to corresponding output file
4. Displays results on console

### Input Format

**Basic Format** (input1.txt):
```
n                           # Number of data points
x₀ y₀                      # First data point
x₁ y₁                      # Second data point
...
xₙ₋₁ yₙ₋₁                  # Last data point
m                           # Number of interpolation points
x_interp₁                   # First interpolation point
x_interp₂                   # Second interpolation point
...
x_interpₘ                   # Last interpolation point
```

**With Error Analysis** (input2.txt):
```
n                           # Number of initial data points
x₀ y₀
x₁ y₁
...
xₙ₋₁ yₙ₋₁
m                           # Number of interpolation points
x_interp₁
... 
x_interpₘ
xₙ yₙ                      # Additional point for error comparison
```

**Example** (Forward/Backward Interpolation):
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

**Example** (Divided Difference with irregular spacing):
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

### Output Format

Each run produces:
1. **Header** with method name
2. **Data points table** showing input data
3. **Difference table** (forward/backward/divided)
4. **Interpolation results** for each requested point
5. **Optional error analysis** if additional point provided

### Customization Options

**Toggle Intermediate Output**:
```cpp
bool printIntermediate = true;  // Set to false for final results only
```

**Adjust Precision**:
```cpp
fout << fixed << setprecision(6);  // Change 6 to desired decimal places
```

**Modify Tolerance**:
```cpp
const double EPSILON = 1e-12;  // Adjust for your numerical requirements
```

---

## 📚 References

- **Numerical Methods for Engineers** by Chapra & Canale
- [Wikipedia: Polynomial Interpolation](https://en.wikipedia.org/wiki/Polynomial_interpolation)

---

## 👨‍💻 Author

**Abir Hasan Arko**  
Roll: 2207053   
CSE, KUET    
[![GitHub](https://img.shields.io/badge/GitHub-AbirHasanArko-181717?style=flat&logo=github)](https://github.com/AbirHasanArko)

---

<div align="center">

**[⬆ Back to Top](#interpolation-and-approximation)**

Part of the [Numerical Computing Suite](https://github.com/AbirHasanArko/Numerical-Computing-Suite)

</div>