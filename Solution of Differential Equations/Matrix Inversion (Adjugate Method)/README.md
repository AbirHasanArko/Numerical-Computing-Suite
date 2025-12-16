# Matrix Inversion (Adjugate Method)

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](matrix-inversion-adjugate.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

---

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Theory](#-theory)
  - [Adjugate / Adjoint Matrix](#adjugate--adjoint-matrix)
  - [Inverse Formula](#inverse-formula)
  - [When Inverse Exists](#when-inverse-exists)
- [Algorithm](#-algorithm)
- [Complexity](#-complexity)
- [Input & Output Format](#-input--output-format)
- [Compilation](#-compilation)
- [References](#-references)

---

## 📖 Introduction

The **Adjugate (Adjoint) Method** computes the inverse of a square matrix using:

- determinant, **det(A)**
- cofactor matrix
- adjugate matrix **adj(A)** (transpose of cofactor matrix)

This is a popular **teaching method** because it matches the theory taught in class.  
However, it becomes slow for large matrices, so it is mainly used for **small n**.

---

## 🧠 Theory

### Adjugate / Adjoint Matrix

For a square matrix **A**, the **cofactor** of element \(a_{ij}\) is:

\[
C_{ij} = (-1)^{i+j} M_{ij}
\]

where \(M_{ij}\) is the **minor determinant** (determinant of the matrix after removing row i and column j).

The **cofactor matrix** is \([C_{ij}]\).

The **adjugate** matrix is:

\[
\text{adj}(A) = (\text{cofactor}(A))^T
\]

### Inverse Formula

If \(\det(A) \neq 0\),

\[
A^{-1} = \frac{\text{adj}(A)}{\det(A)}
\]

### When Inverse Exists

- If **det(A) = 0** → matrix is **singular** → inverse does **not** exist  
- If **det(A) ≠ 0** → inverse exists

---

## 🧩 Algorithm

1. Read \(n\) and the \(n \times n\) matrix \(A\)
2. Compute \(\det(A)\)
3. If \(|\det(A)|\) is near zero, report **No Inverse**
4. Compute all cofactors \(C_{ij}\)
5. Form adjugate: \(\text{adj}(A) = C^T\)
6. Compute inverse: \(A^{-1} = \frac{\text{adj}(A)}{\det(A)}\)
7. Print: input matrix, det(A), cofactor matrix, adjugate, inverse

---

## ⏱️ Complexity

This implementation uses **cofactor expansion** for determinant/minors:

- Determinant (Laplace expansion): grows very fast (≈ **O(n!)**)
- Works best for **n ≤ 6** (recommended for lab)

For real engineering problems, methods like **Gauss-Jordan** or **LU** are preferred.

---

## 🔢 Input & Output Format

### Input (`input.txt`)
```
n
a11 a12 ... a1n
a21 a22 ... a2n
...
an1 an2 ... ann
```

### Output (`output.txt`)
- Prints the matrix
- det(A)
- cofactor matrix
- adjugate matrix
- inverse matrix (if exists)

---

## 🛠️ Compilation

```bash
g++ -std=c++17 -O2 matrix-inversion-adjugate.cpp -o inv
./inv
```

---

## 📚 References

- Kreyszig — *Advanced Engineering Mathematics*
- Burden & Faires — *Numerical Analysis*
- Wikipedia: Adjugate matrix, Cofactor expansion
