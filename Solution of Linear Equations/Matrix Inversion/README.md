# Matrix Inversion (Adjugate Method)

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](matrix-inversion.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

---

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Theory & Algorithm](#-theory--algorithm)
  - [Mathematical Foundation](#mathematical-foundation)
  - [Algorithm Steps](#algorithm-steps)
  - [Singularity Detection](#singularity-detection)
  - [Complexity Analysis](#complexity-analysis)
- [Implementation Details](#-implementation-details)
  - [Key Components](#key-components)
  - [Determinant & Adjoint Calculation](#determinant--adjoint-calculation)
  - [Numerical Stability](#numerical-stability)
- [Complete C++ Implementation](#-complete-c-implementation)
- [Usage Examples](#-usage-examples)
  - [Example: 3x3 Matrix](#example-3x3-matrix)
- [Compilation and Execution](#-compilation-and-execution)
- [Applications](#-applications)
- [References](#-references)

---

## 📖 Introduction

The **Matrix Inversion by Adjugate Method** is a classical approach to finding the inverse of a square matrix. It leverages the properties of determinants and cofactors to compute the inverse, provided the matrix is non-singular (i.e., has a non-zero determinant).

This C++ implementation reads a square matrix from a file, computes its inverse using the adjugate (adjoint) and determinant, and writes the result to an output file. If the matrix is singular, it reports accordingly.

### Features
- ✅ Handles any square matrix (n x n)
- ✅ Detects singular (non-invertible) matrices
- ✅ High precision output (configurable)
- ✅ File-based I/O for easy testing
- ✅ Modular code (determinant, adjoint, inverse)

---

## 🧮 Theory & Algorithm

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

## 💻 Implementation Details

### Key Components

- **Determinant Function**: Recursively computes determinant via cofactor expansion.
- **Adjoint Function**: Computes the matrix of cofactors and transposes it.
- **Inverse Function**: Combines determinant and adjoint to compute the inverse.
- **File I/O**: Reads matrix from `input.txt`, writes result to `output.txt`.
- **Precision**: Uses `setprecision(6)` for output.

### Determinant & Adjoint Calculation

- **Determinant**: For $n=1$ or $n=2$, uses direct formula. For $n>2$, expands along the first row.
- **Adjoint**: For each element, computes its cofactor (minor with sign), then transposes.

### Numerical Stability

- Uses a small threshold ($\varepsilon = 10^{-12}$) to detect near-zero determinants.
- All calculations use `double` for higher precision.

---

## 🔧 Complete C++ Implementation

See [matrix-inversion.cpp](matrix-inversion.cpp) for the full code.

---

## 📊 Usage Examples

### Input Format

```
n
A₁₁ A₁₂ ... A₁ₙ
A₂₁ A₂₂ ... A₂ₙ
...
Aₙ₁ Aₙ₂ ... Aₙₙ
```

### Example: 3x3 Matrix

**Input (input.txt):**
```
3
4 7 2
3 6 1
2 5 3
```

**Output (output.txt):**
```
Inverse using adj(A)/det(A):
    1.444444   -1.222222   -0.555556
   -0.777778    0.888889    0.222222
    0.333333   -0.666667    0.333333
```

**If the matrix is singular:**
```
Matrix is singular.
```

---

## 🎯 Compilation and Execution

### Compile
```bash
g++ -std=c++17 -O2 matrix-inversion.cpp -o matrix-inv
```

### Run
```bash
./matrix-inv
```

### Requirements
- C++11 or later
- Input file: `input.txt` in the same directory
- Output file: `output.txt` (automatically created)

---

## 🔬 Applications

Matrix inversion is fundamental in:

- **Solving systems of linear equations**
- **Computer graphics** (transformations)
- **Control systems** (state-space analysis)
- **Cryptography** (Hill cipher)
- **Statistics** (covariance matrix inversion)
- **Engineering simulations**

---

## 📚 References

- [Matrix Inversion - Wikipedia](https://en.wikipedia.org/wiki/Invertible_matrix)
- Numerical Methods For Engineers by Raymond Canale and Steven C. Chapra

---

## 👤 Author

**Part of the [Numerical Computing Suite](../../) by [AbirHasanArko](https://github.com/AbirHasanArko)**   
Roll: 2207053   
Department of CSE, KUET   
