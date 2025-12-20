# Solution of Non-linear Equations

[![C++](https://img.shields.io/badge/Language-C++-00599C?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)
[![Numerical Methods](https://img.shields.io/badge/Topic-Numerical%20Methods-FF6B6B?style=for-the-badge)](https://github.com/AbirHasanArko/Numerical-Computing-Suite)

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Overview of Non-linear Equation Methods](#-overview-of-non-linear-equation-methods)
- [Mathematical Foundation](#-mathematical-foundation)
  - [What is a Non-linear Equation?](#what-is-a-non-linear-equation)
  - [Root-finding Problem](#root-finding-problem)
  - [Convergence and Error](#convergence-and-error)
- [Methods Implemented](#-methods-implemented)
  - [1. Bisection Method](#1-bisection-method)
  - [2. False-Position (Regula Falsi) Method](#2-false-position-regula-falsi-method)
  - [3. Newton-Raphson Method](#3-newton-raphson-method)
  - [4. Secant Method](#4-secant-method)
- [Method Comparison](#-method-comparison)
- [Applications](#-applications)
- [Implementation Structure](#-implementation-structure)
- [Getting Started](#-getting-started)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

This section of the **Numerical Computing Suite** focuses on **numerical methods for solving non-linear equations**—that is, finding roots of equations where the function is not a simple straight line or polynomial of degree one. Such equations arise in almost every field of science and engineering, and analytical solutions are often impossible. Numerical methods provide robust, iterative approaches to approximate these roots to any desired accuracy.

**The Core Problem**: Given a function $f(x)$, find $x^*$ such that $f(x^*) = 0$.

Root-finding is fundamental in:
- **Physics**: Equilibrium points, energy levels
- **Engineering**: Circuit analysis, stress calculations
- **Mathematics**: Transcendental equations, optimization
- **Finance**: Internal rate of return, option pricing

This collection provides four classic and widely used methods, each with unique strengths and convergence properties.

---

## 🔍 Overview of Non-linear Equation Methods

This repository implements four major root-finding algorithms:

| Method         | Bracketing? | Derivative Needed? | Order of Convergence | Robustness      |
|---------------|-------------|--------------------|---------------------|-----------------|
| **Bisection** | Yes         | No                 | Linear              | Very robust     |
| **False-Position** | Yes    | No                 | Linear (faster than bisection) | Robust      |
| **Newton-Raphson** | No     | Yes                | Quadratic           | Fast, less robust |
| **Secant**    | No          | No                 | Superlinear         | Fast, less robust |

---

## 📐 Mathematical Foundation

### What is a Non-linear Equation?

A **non-linear equation** is any equation where the variable $x$ appears with exponents other than one, inside transcendental functions (sin, exp, log, etc.), or in products with itself. Examples:
- $x^3 - x - 2 = 0$
- $e^x + x = 5$
- $\sin(x) = x/2$

### Root-finding Problem

The goal is to find $x^*$ such that $f(x^*) = 0$. Analytical solutions are rare, so we use iterative numerical methods:
- **Bracketing methods**: Start with an interval $[a, b]$ where $f(a)$ and $f(b)$ have opposite signs (guaranteed root by Intermediate Value Theorem).
- **Open methods**: Use one or two initial guesses and iterate using function values (and possibly derivatives).

### Convergence and Error

- **Convergence**: How quickly the method approaches the root.
- **Order of convergence**: Linear (slow), superlinear, quadratic (fast).
- **Stopping criteria**: $|f(x_n)| < \epsilon$ or $|x_{n+1} - x_n| < \epsilon$ for some small $\epsilon$.
- **Error estimation**: Each method provides a way to estimate or bound the error at each step.

---

## 🛠️ Methods Implemented

### 1. Bisection Method

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-blue?style=for-the-badge)](./Bisection%20Method/)

- **Bracketing method**: Requires $f(a) \cdot f(b) < 0$.
- Repeatedly halves the interval $[a, b]$.
- Always converges, but slowly (linear rate).
- Very robust, guaranteed to find a root if one exists in $[a, b]$.

### 2. False-Position (Regula Falsi) Method

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-green?style=for-the-badge)](./False-Position%20Method/)

- **Bracketing method**: Like bisection, but uses a secant line to estimate the root.
- Faster than bisection in many cases.
- Still robust, but can stagnate if one endpoint doesn't move.

### 3. Newton-Raphson Method

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-orange?style=for-the-badge)](./Newton-Raphson%20Method/)

- **Open method**: Needs a good initial guess and $f'(x)$.
- Iterates: $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$
- Quadratic convergence (very fast) if initial guess is close.
- Can fail if $f'(x)$ is zero or guess is poor.

### 4. Secant Method

[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-red?style=for-the-badge)](./Secant%20Method/)

- **Open method**: Uses two initial guesses, no derivative needed.
- Iterates using a secant line through last two points.
- Superlinear convergence (faster than bisection, slower than Newton-Raphson).
- Can fail if guesses are poor or function is ill-behaved.

---

## 📊 Method Comparison

| Aspect              | Bisection | False-Position | Newton-Raphson | Secant |
|---------------------|-----------|----------------|----------------|--------|
| **Bracketing**      | Yes       | Yes            | No             | No     |
| **Derivative**      | No        | No             | Yes            | No     |
| **Convergence**     | Linear    | Linear         | Quadratic      | Superlinear |
| **Robustness**      | High      | High           | Medium         | Medium |
| **Speed**           | Slow      | Moderate       | Fast           | Fast   |
| **Failure Modes**   | None (if root exists) | Stagnation | Divergence, division by zero | Divergence |

---

## 🎯 Applications

- **Engineering**: Solving circuit equations, stress analysis, fluid flow
- **Physics**: Finding equilibrium, energy levels, resonance frequencies
- **Mathematics**: Solving transcendental and polynomial equations
- **Finance**: Calculating IRR, bond yields, option pricing
- **Computer Science**: Non-linear optimization, graphics, machine learning

---

## 📁 Implementation Structure

Each method folder contains:

```
Method Name/
├── README.md                      # Detailed theory and examples
├── [method-filename].cpp          # C++ implementation
├── input.txt                      # Test case(s)
├── output.txt                     # Expected output
```

### Common Features Across All Implementations

- ✅ File-based I/O for reproducibility
- ✅ Multiple test cases per file
- ✅ Formatted output with configurable precision
- ✅ Step-by-step iteration logs
- ✅ Error and convergence checks
- ✅ Input validation and meaningful error messages
- ✅ Well-commented, modular code

---

## 🚀 Getting Started

### Prerequisites
- **C++ Compiler**: g++, clang++, or MSVC
- **C++ Standard**: C++11 or later
- **Text Editor/IDE**: VS Code, CLion, or any preferred editor
- **Basic Knowledge**: Calculus, function analysis

### Compilation

Navigate to any method folder and compile:

```bash
# Standard compilation
g++ -o solve [method-filename].cpp -std=c++11

# With optimization (recommended)
g++ -o solve [method-filename].cpp -std=c++17 -O2

# With warnings (for development)
g++ -o solve [method-filename].cpp -std=c++17 -O2 -Wall -Wextra
```

### Execution

```bash
# Run the compiled program
./solve

# On Windows
solve.exe
```

The program:
1. Reads input from `input.txt`
2. Processes all root-finding requests
3. Writes results to `output.txt`
4. Displays results on console

### Input Format

**Basic Format** (input.txt):
```
function_expression         # e.g., x^3 - x - 2
interval_or_guesses         # e.g., 1 2 (for bracketing) or 1 2 (for secant)
max_iterations              # e.g., 100
error_tolerance             # e.g., 1e-6
```

**Example** (Bisection/False-Position):
```
x^3 - x - 2
1 2
100
1e-6
```

**Example** (Newton-Raphson):
```
x^3 - x - 2
1.5
100
1e-6
```

### Output Format

Each run produces:
1. **Header** with method name
2. **Iteration table** showing $x_n$, $f(x_n)$, error
3. **Final root estimate** and number of iterations
4. **Convergence status**

### Customization Options

- **Adjust Precision**:
  ```cpp
  fout << fixed << setprecision(8);  // Change 8 to desired decimal places
  ```
- **Modify Tolerance**:
  ```cpp
  const double EPSILON = 1e-8;  // Adjust for your numerical requirements
  ```

---

## 📚 References

- **Numerical Methods for Engineers** by Chapra & Canale
- [Wikipedia: Root-finding Algorithms](https://en.wikipedia.org/wiki/Root-finding_algorithms)

---

## 👨‍💻 Author

**MD. Shomik Shahriar**  
Roll: 2207041   
CSE, KUET    
[![GitHub](https://img.shields.io/badge/GitHub-HapiGuy-181717?style=flat&logo=github)](https://github.com/Hapi-Guy)

---

<div align="center">

**[⬆ Back to Top](#solution-of-non-linear-equations)**

Part of the [Numerical Computing Suite](https://github.com/AbirHasanArko/Numerical-Computing-Suite)

</div>
