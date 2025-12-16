# Curve Fitting (Regression)

[![C++](https://img.shields.io/badge/Language-C++-00599C?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)
[![Numerical Methods](https://img.shields.io/badge/Topic-Numerical%20Methods-FF6B6B?style=for-the-badge)](https://en.wikipedia.org/wiki/Curve_fitting)

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Overview](#-overview)
- [Methods Included](#-methods-included)
  - [1. Least Squares Line (Linear Regression)](#1-least-squares-line-linear-regression)
  - [2. Least Squares Polynomial (Quadratic)](#2-least-squares-polynomial-quadratic)
  - [3. Non-Linear Curve Fitting](#3-non-linear-curve-fitting)
- [Implementation Structure](#-implementation-structure)
- [Getting Started](#-getting-started)
- [References](#-references)

---

## 📖 Introduction

This folder contains C++ implementations commonly used in **Numerical Methods** courses for solving
**curve fitting and regression problems** based on experimental or tabulated data.

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

---

## 🔍 Overview

| Method | Curve Type | Main Idea | Notes |
|---|---|---|---|
| Least Squares Line | Linear | Minimize squared errors | Simple and widely used |
| Least Squares Polynomial | Quadratic | Higher-degree fitting | Better for curved data |
| Non-Linear Curve Fitting | Exponential / Power | Transform to linear form | Model-dependent accuracy |

---

## ✅ Methods Included

### 1. Least Squares Line (Linear Regression)

Fits a straight-line equation of the form:

\[
y = a + bx
\]

The constants \(a\) and \(b\) are obtained by solving the **normal equations** derived from the least squares principle:

\[
\sum y = na + b\sum x
\]

\[
\sum xy = a\sum x + b\sum x^2
\]

This method is suitable when the data shows a **linear relationship** between variables.

---

### 2. Least Squares Polynomial (Quadratic)

Fits a second-degree polynomial curve of the form:

\[
y = a + bx + cx^2
\]

The coefficients \(a\), \(b\), and \(c\) are determined by solving the following system of normal equations:

\[
\sum y = na + b\sum x + c\sum x^2
\]

\[
\sum xy = a\sum x + b\sum x^2 + c\sum x^3
\]

\[
\sum x^2y = a\sum x^2 + b\sum x^3 + c\sum x^4
\]

This method is used when the data follows a **curved trend** that cannot be accurately modeled by a straight line.

---

### 3. Non-Linear Curve Fitting

Non-linear curve fitting is used when the relationship between variables is inherently non-linear.
Common non-linear models include:

- Exponential curve:  
  \[
  y = ae^{bx}
  \]

- Power curve:  
  \[
  y = ax^b
  \]

These equations are transformed into **linear form** using logarithmic transformations, after which
the least squares method is applied to determine the constants.

This method is widely used in **engineering, physics, and biological data analysis**.

---

## Implementation Structure

Curve Fitting/
├─ README.md
├─ Least Squares Line/
│  ├─ least-squares-line.cpp
│  ├─ input.txt
│  └─ output.txt
├─ Least Squares Polynomial/
│  ├─ least-squares-polynomial.cpp
│  ├─ input.txt
│  └─ output.txt
└─ Non-Linear Curve Fitting/
├─ non-linear-curve-fitting.cpp
├─ input.txt
└─ output.txt

---

## Getting Started

1. Navigate to the desired method folder.
2. Check `input.txt` for the required input format.
3. Compile and run the program:

```bash
g++ -std=c++17 -O2 -o main main.cpp
./main
```
---

## 📚 References

- S. S. Sastry, *Introductory Methods of Numerical Analysis*  
- B. S. Grewal, *Numerical Methods in Engineering and Science*
