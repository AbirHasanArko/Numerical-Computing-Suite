# Runge-Kutta Method (RK4)

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](runge-kutta-rk4.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

---

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Theory](#-theory)
  - [Initial Value Problem (IVP)](#initial-value-problem-ivp)
  - [RK4 Formula](#rk4-formula)
- [Algorithm](#-algorithm)
- [Complexity](#-complexity)
- [Input & Output Format](#-input--output-format)
- [How to Change the Function f(x,y)](#-how-to-change-the-function-fxy)
- [Compilation](#-compilation)
- [References](#-references)

---

## 📖 Introduction

The **Runge–Kutta 4th Order method (RK4)** is one of the most widely used numerical techniques to solve first-order ordinary differential equations of the form:

\[
\frac{dy}{dx} = f(x,y), \quad y(x_0)=y_0
\]

RK4 improves accuracy by taking a weighted average of slopes inside each step.

---

## 🧠 Theory

### Initial Value Problem (IVP)

Given:
- ODE: \(y' = f(x,y)\)
- Initial condition: \(y(x_0)=y_0\)
- Step size: \(h\)

We compute \(y\) at \(x_0+h, x_0+2h, \dots\) up to the target point.

### RK4 Formula

For each step:

\[
k_1 = h f(x_n, y_n)
\]
\[
k_2 = h f(x_n + \tfrac{h}{2}, y_n + \tfrac{k_1}{2})
\]
\[
k_3 = h f(x_n + \tfrac{h}{2}, y_n + \tfrac{k_2}{2})
\]
\[
k_4 = h f(x_n + h, y_n + k_3)
\]

Update:

\[
y_{n+1} = y_n + \frac{k_1 + 2k_2 + 2k_3 + k_4}{6}
\]
\[
x_{n+1} = x_n + h
\]

---

## 🧩 Algorithm

1. Read \(x_0, y_0\)
2. Read target \(x_n\) and step size \(h\)
3. Repeat while \(x < x_n\):
   - compute \(k_1, k_2, k_3, k_4\)
   - update \(y\) and \(x\)
4. Print a table of steps and the final answer

---

## ⏱️ Complexity

If total steps = \(N\), then:

- Time: **O(N)**
- Space: **O(1)** (only a few variables)

---

## 🔢 Input & Output Format

### Input (`input.txt`)
This program provides **built-in examples** for \(f(x,y)\).  
You choose a function by entering an option number.

```
option
x0 y0
x_target h
```

### Functions included
1. \(f(x,y) = x + y\)
2. \(f(x,y) = x - y\)
3. \(f(x,y) = y - x^2 + 1\)
4. \(f(x,y) = x*y\)

You can add more very easily.

---

## ✍️ How to change the function f(x,y)

Open `runge-kutta-rk4.cpp` and edit the function:

```cpp
double f(int option, double x, double y) { ... }
```

Add your own case and select it from input.

---

## 🛠️ Compilation

```bash
g++ -std=c++17 -O2 runge-kutta-rk4.cpp -o rk4
./rk4
```

---

## 📚 References

- Burden & Faires — *Numerical Analysis*
- Kreyszig — *Advanced Engineering Mathematics*
- Wikipedia: Runge–Kutta methods
