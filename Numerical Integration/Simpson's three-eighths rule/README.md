
# Simpson's Three-Eighths (3/8) Rule

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](simpsons-three-eighths-rule.cpp)
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
Simpson’s 3/8 Rule approximates the definite integral by fitting **cubic (3rd-degree polynomials)** over groups of **three consecutive subintervals**.
It provides good accuracy for smooth functions when the number of subintervals is a multiple of three.

---

## 📌 Mathematical Formula
Divide the interval \([a,b]\) into **n** equal subintervals (n must be a multiple of 3).

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
\int_a^b f(x)\,dx \approx \frac{3h}{8}\Big[
y_0 + y_n
+ 3(y_1 + y_2 + y_4 + y_5 + \dots)
+ 2(y_3 + y_6 + y_9 + \dots)
\Big]
\]

where \(y_i = f(x_i)\).

---

## ✅ Validity Condition
- `n` must be a **multiple of 3**
- data must be **equally spaced**

---

## 🧾 Algorithm Steps
1. Read `n`
2. Read `a` and `b`
3. Read `y0..yn` (total `n+1` values)
4. Check `n` is divisible by 3 (otherwise print error for that case)
5. Compute `h = (b-a)/n`
6. Compute:
   - `sum_3 = y1 + y2 + y4 + y5 + ...` (indices not divisible by 3)
   - `sum_2 = y3 + y6 + y9 + ...` (indices divisible by 3, excluding endpoints)
7. Compute integral using Simpson’s 3/8 formula
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
