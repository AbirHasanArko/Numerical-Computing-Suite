# Least Squares Line (Linear Regression)

---

## 📖 Introduction
The least squares line is used to find the **best fitting straight line**
for a given set of experimental data points.
The equation of the straight line is assumed to be:

\[
y = a + bx
\]

The constants `a` and `b` are determined such that the sum of squares of
errors between observed and computed values is minimum.

---

## 📌 Mathematical Formula
Normal equations for linear regression:

\[
\sum y = na + b\sum x
\]

\[
\sum xy = a\sum x + b\sum x^2
\]

Solving these equations gives the values of `a` and `b`.

---

## 🧾 Algorithm Steps
1. Read the number of data points `n`
2. Read values of `x` and `y`
3. Compute required summations:
   - Σx, Σy, Σxy, Σx²
4. Solve the normal equations
5. Obtain constants `a` and `b`
6. Print the fitted equation

---

## ⚙️ Implementation Notes
- Applicable only for linear relationships
- Simple and widely used regression method

---

## 🧪 Usage Example
Refer to `input.txt` and `output.txt`.

---

## 📚 References
- S. S. Sastry  
- B. S. Grewal