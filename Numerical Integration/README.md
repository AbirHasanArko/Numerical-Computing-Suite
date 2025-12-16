# Numerical Integration (Simpson's Rules)

[![C++](https://img.shields.io/badge/Language-C++-00599C?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)
[![Numerical Methods](https://img.shields.io/badge/Topic-Numerical%20Integration-2F855A?style=for-the-badge)](#)
[![Simpson Rules](https://img.shields.io/badge/Method-Simpson's%20Rules-6B46C1?style=for-the-badge)](#)

## 📑 Table of Contents
- [📖 Introduction](#-introduction)
- [🧠 Methods Included](#-methods-included)
- [📌 Input Format](#-input-format)
- [⚙️ How to Compile & Run](#-how-to-compile--run)
- [✅ Notes & Validity Conditions](#-notes--validity-conditions)
- [📚 References](#-references)

---

## 📖 Introduction
**Numerical Integration** (Numerical Quadrature) is used to approximate definite integrals when:
- an analytical integral is difficult or impossible,
- the function values are given in **tabular form**,
- or a quick approximation is required.

Simpson’s rules are popular Newton–Cotes methods that assume **equally spaced** x-values and approximate the integrand using low-degree polynomials.

---

## 🧠 Methods Included
- **Simpson's One-Third (1/3) Rule**
- **Simpson's Three-Eighths (3/8) Rule**

Each method folder contains:
- `README.md` (theory + algorithm + notes)
- C++ source code (`.cpp`)
- `input.txt` (multiple sample test cases)
- `output.txt` (sample output format)

---

## 📌 Input Format
Both programs work with **tabulated y-values** (not a hardcoded function).  
Each test case is given as:

```
n
a b
y0 y1 y2 ... yn
```

Where:
- `n` = number of subintervals
- `a, b` = lower and upper limits
- `y0..yn` = function values at equally spaced points

Multiple test cases can be placed one after another in the same input file (the program reads until EOF).

---

## ⚙️ How to Compile & Run
Example (Linux / Mac):
```bash
g++ simpsons-one-third-rule.cpp -o s13
./s13 < input.txt
```

Windows (MinGW):
```bash
g++ simpsons-one-third-rule.cpp -o s13.exe
s13.exe < input.txt
```

---

## ✅ Notes & Validity Conditions
- The x-values must be **equally spaced**, i.e., constant step size `h = (b-a)/n`.
- **Simpson’s 1/3 Rule** requires `n` to be **even**.
- **Simpson’s 3/8 Rule** requires `n` to be a **multiple of 3**.

---

## 📚 References
- S. S. Sastry, *Introductory Methods of Numerical Analysis*
- B. S. Grewal, *Numerical Methods in Engineering and Science*
- Any standard Numerical Methods course notes

---

<div align="center">

**[⬆ Back to Top](#numerical-integration-simpsons-rules)**

</div>
