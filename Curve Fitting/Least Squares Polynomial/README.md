# Least Squares Polynomial (Quadratic Curve)

---

## 📖 Introduction
Polynomial curve fitting is used when the relationship between variables
is **non-linear**.
A second-degree polynomial is assumed in the form:

\[
y = a + bx + cx^2
\]

The coefficients are determined using the least squares method.

---

## 📌 Mathematical Formula
Normal equations for quadratic fitting:

\[
\sum y = na + b\sum x + c\sum x^2
\]

\[
\sum xy = a\sum x + b\sum x^2 + c\sum x^3
\]

\[
\sum x^2y = a\sum x^2 + b\sum x^3 + c\sum x^4
\]

Solving these equations gives `a`, `b`, and `c`.

---

## 🧾 Algorithm Steps
1. Read number of data points
2. Read `x` and `y` values
3. Compute required summations
4. Form the system of normal equations
5. Solve the equations
6. Display polynomial coefficients

---

## ⚙️ Implementation Notes
- Suitable for curved data trends
- Higher accuracy than linear regression for non-linear data

---

## 🧪 Usage Example
Refer to `input.txt` and `output.txt`.

---

## 📚 References
- S. S. Sastry  
- B. S. Grewal