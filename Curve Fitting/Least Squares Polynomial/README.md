## Least Squares Polynomial (Quadratic Curve)

## 📑 Table of Contents
- [Introduction](#-introduction-1)
- [Mathematical Model](#-mathematical-model-1)
- [Normal Equations](#-normal-equations-1)
- [Algorithm](#-algorithm-1)
- [Applications](#-applications-1)

---

## 📖 Introduction
When experimental data does not follow a straight-line pattern,
a polynomial curve provides a better approximation.
In this method, a **second-degree polynomial** is fitted using the least squares principle.

---

## 📌 Mathematical Model
The quadratic polynomial is assumed as:

\[
y = a + bx + cx^2
\]

where \(a\), \(b\), and \(c\) are constants.

---

## 📐 Normal Equations
The normal equations are:

\[
\sum y = na + b\sum x + c\sum x^2
\]

\[
\sum xy = a\sum x + b\sum x^2 + c\sum x^3
\]

\[
\sum x^2y = a\sum x^2 + b\sum x^3 + c\sum x^4
\]

Solving these equations gives the polynomial coefficients.

---

## 🧾 Algorithm
1. Read the number of observations
2. Read the values of \(x\) and \(y\)
3. Compute required summations
4. Form the normal equations
5. Solve for \(a\), \(b\), and \(c\)
6. Display the fitted polynomial equation

---

## 🧪 Applications
- Non-linear experimental data
- Engineering curve modeling
- Scientific data approximation