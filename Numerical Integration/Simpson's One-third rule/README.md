# Simpson's One-Third (1/3) Rule

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](simpsons-one-third-rule.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

---

## 📑 Table of Contents
- [📖 Introduction](#-introduction)
- [📌 Mathematical Formula](#-mathematical-formula)
- [✅ Validity Condition](#-validity-condition)
- [🧾 Algorithm Steps](#-algorithm-steps)
- [⚙️ Implementation Notes](#-implementation-notes)
- [🧪 Usage Example](#-usage-example)
- [📚 References](#-references)

---

## 📖 Introduction
Simpson’s 1/3 Rule approximates the integral by fitting **parabolas (2nd-degree polynomials)** over pairs of subintervals.
It is generally more accurate than the trapezoidal rule for smooth functions.

---

## 📌 Mathematical Formula
Divide the interval \([a,b]\) into **n** equal subintervals (n must be even).

- Step size:  
\(
h = \frac{b-a}{n}
\)

- Points:  
\(
x_i = a + ih,\;\; i = 0,1,\dots,n
\)

Then:

\[
\int_a^b f(x)\,dx \approx \frac{h}{3}\Big[
y_0 + y_n + 4(y_1+y_3+\dots+y_{n-1}) + 2(y_2+y_4+\dots+y_{n-2})
\Big]
\]

where \(y_i = f(x_i)\).

---

## ✅ Validity Condition
- `n` must be **even**
- data must be **equally spaced**

---

## 🧾 Algorithm Steps
1. Read `n`
2. Read `a` and `b`
3. Read `y0..yn` (total `n+1` values)
4. Check `n` is even (otherwise print error for that case)
5. Compute `h = (b-a)/n`
6. Compute:
   - `sum_odd = y1 + y3 + ... + y(n-1)`
   - `sum_even = y2 + y4 + ... + y(n-2)`
7. Compute integral using Simpson’s 1/3 formula
8. Print result in a clear formatted way

---

## ⚙️ Implementation Notes
- The program supports **multiple test cases** until EOF.
- It prints the weights and partial sums to match typical lab-output style.

---

## 🧪 Usage Example
See `input.txt` and `output.txt` inside this folder.

---

## 📚 References
- S. S. Sastry, *Introductory Methods of Numerical Analysis*
- B. S. Grewal, *Numerical Methods in Engineering and Science*
