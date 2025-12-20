# Runge-Kutta Method (RK4)

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](runge-kutta-rk4.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

---

## 📑 Table of Contents

- [📖 Introduction](#-introduction)
- [🧠 Theory](#-theory)
  - [Initial Value Problem (IVP)](#initial-value-problem-ivp)
  - [RK4 Formula](#rk4-formula)
- [🧩 Algorithm](#-algorithm)
- [⏱️ Complexity](#-complexity)
- [🔢 Input & Output Format](#-input--output-format)
- [✍️ How to Change the Function f(x,y)](#-how-to-change-the-function-fxy)
- [🛠️ Compilation](#-compilation)
- [📚 References](#-references)

---

## 📖 Introduction

The **Runge–Kutta 4th Order method (RK4)** is one of the most widely used numerical techniques to solve first-order ordinary differential equations of the form:

```
dy/dx = f(x, y),    y(x₀) = y₀
```

RK4 improves accuracy by taking a weighted average of slopes inside each step.

---

## 🧠 Theory

### Initial Value Problem (IVP)

Given:
- ODE: `y' = f(x, y)`
- Initial condition: `y(x₀) = y₀`
- Step size: `h`

We compute `y` at `x₀ + h`, `x₀ + 2h`, … up to the target point.

### RK4 Formula

For each step:

```
k₁ = h * f(xₙ,      yₙ)
k₂ = h * f(xₙ+h/2,  yₙ+k₁/2)
k₃ = h * f(xₙ+h/2,  yₙ+k₂/2)
k₄ = h * f(xₙ+h,    yₙ+k₃)
```

Update:

```
yₙ₊₁ = yₙ + (k₁ + 2k₂ + 2k₃ + k₄) / 6
xₙ₊₁ = xₙ + h
```

---

## 🧩 Algorithm

1. Read `x₀`, `y₀`
2. Read target `xₙ` and step size `h`
3. Repeat while `x < xₙ`:
   - compute `k₁`, `k₂`, `k₃`, `k₄`
   - update `y` and `x`
4. Print a table of steps and the final answer

---

## ⏱️ Complexity

If total steps = `N`:

- Time: **O(N)**
- Space: **O(1)** (only a few variables needed)

---

## 🔢 Input & Output Format

### Input (`input.txt`)
The program provides **built-in examples** for `f(x, y)`.  
You choose a function by entering an option number.

```
option
x0 y0
x_target h
```

### Functions included
1. `f(x, y) = x + y`
2. `f(x, y) = x - y`
3. `f(x, y) = y - x^2 + 1`
4. `f(x, y) = x * y`

You can add more functions easily.

---

## ✍️ How to Change the Function f(x, y)

Open `runge-kutta-rk4.cpp` and edit the function:

```cpp
double f(int option, double x, double y) { ... }
```

Add your own `case` and select it from the input.

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