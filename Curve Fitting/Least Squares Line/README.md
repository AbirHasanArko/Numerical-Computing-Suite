## Least Squares Line (Linear Regression)

## 📑 Table of Contents
- [Introduction](#-introduction)
- [Mathematical Model](#-mathematical-model)
- [Normal Equations](#-normal-equations)
- [Algorithm](#-algorithm)
- [Applications](#-applications)

---

## 📖 Introduction
The Least Squares Line method is used to find the **best fitting straight line**
for a given set of experimental or observed data points.
The method minimizes the **sum of squares of vertical errors**
between the observed values and the computed values.

---

## 📌 Mathematical Model
The straight line equation is assumed as:

\[
y = a + bx
\]

where \(a\) and \(b\) are constants to be determined.

---

## 📐 Normal Equations
Using the least squares principle, the normal equations are:

\[
\sum y = na + b\sum x
\]

\[
\sum xy = a\sum x + b\sum x^2
\]

Solving these equations gives the values of \(a\) and \(b\).

---

## 🧾 Algorithm
1. Read the number of data points \(n\)
2. Read values of \(x\) and \(y\)
3. Compute \(\sum x\), \(\sum y\), \(\sum xy\), \(\sum x^2\)
4. Form the normal equations
5. Solve for \(a\) and \(b\)
6. Display the fitted line equation

---

## 🧪 Applications
- Linear trend analysis
- Engineering measurements
- Economics and statistical data analysis