# Solution of Differential Equations

[![C++](https://img.shields.io/badge/Language-C++-00599C?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)
[![Numerical Methods](https://img.shields.io/badge/Topic-Numerical%20Methods-FF6B6B?style=for-the-badge)](https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations)

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Overview](#-overview)
- [Methods Included](#-methods-included)
  - [1. Matrix Inversion (Adjugate Method)](#1-matrix-inversion-adjugate-method)
  - [2. Runge-Kutta Method (RK4)](#2-runge-kutta-method-rk4)
- [Implementation Structure](#-implementation-structure)
- [Getting Started](#-getting-started)
- [References](#-references)

---

## 📖 Introduction

This folder contains C++ implementations commonly used in **Numerical Methods** courses for solving problems related to **differential equations** and the **linear algebra** that often appears alongside them (e.g., boundary value problems, coupled systems, discretization).

You will find:
- A **Matrix Inversion** program using the **Adjugate (Adjoint) Matrix** approach.
- A classic **Runge–Kutta 4th Order (RK4)** solver for first-order ODEs.

Each method folder includes:
- `README.md` (theory + algorithm + usage)
- `*.cpp` source file
- `input.txt` and `output.txt` samples

---

## 🔍 Overview

| Method | Problem Type | Main Idea | Notes |
|---|---|---|---|
| Matrix Inversion (Adjugate) | Linear Algebra | \(A^{-1} = \frac{\text{adj}(A)}{\det(A)}\) | Best for small matrices (teaching) |
| Runge-Kutta (RK4) | ODE (Initial Value) | Weighted slope averaging | Accurate and widely used |

---

## ✅ Methods Included

### 1. Matrix Inversion (Adjugate Method)

📁 Folder: **Matrix Inversion (Adjugate Method)**  
Computes determinant, cofactor matrix, adjugate, and finally the inverse of a square matrix (if invertible).

### 2. Runge-Kutta Method (RK4)

📁 Folder: **Runge-Kutta Method (RK4)**  
Solves \(\frac{dy}{dx}=f(x,y)\) with initial condition \(y(x_0)=y_0\) over a step size \(h\).

---

## 🗂️ Implementation Structure

```
Solution of Differential Equations/
├─ README.md
├─ Matrix Inversion (Adjugate Method)/
│  ├─ README.md
│  ├─ matrix-inversion-adjugate.cpp
│  ├─ input.txt
│  └─ output.txt
└─ Runge-Kutta Method (RK4)/
   ├─ README.md
   ├─ runge-kutta-rk4.cpp
   ├─ input.txt
   └─ output.txt
```

---

## 🚀 Getting Started

1. Open a method folder.
2. Check `input.txt` for sample input format.
3. Compile and run:

```bash
g++ -std=c++17 -O2 -o main main.cpp
./main
```

(Or compile using the provided `.cpp` file name.)

---

## 📚 References

- Burden & Faires — *Numerical Analysis*
- Kreyszig — *Advanced Engineering Mathematics*
- Wikipedia: Numerical methods for ODEs, Runge–Kutta methods, Adjugate matrix
