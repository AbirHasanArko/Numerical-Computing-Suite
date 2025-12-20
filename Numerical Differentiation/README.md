# Numerical Differentiation

[![C++](https://img.shields.io/badge/Language-C++-00599C?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)
[![Numerical Methods](https://img.shields.io/badge/Topic-Numerical%20Methods-FF6B6B?style=for-the-badge)](https://github.com/AbirHasanArko/Numerical-Computing-Suite)

## 📑 Table of Contents

- [Introduction](#-introduction)
- [Overview of Numerical Differentiation](#-overview-of-numerical-differentiation)
- [Mathematical Foundation](#-mathematical-foundation)
  - [What is Numerical Differentiation?](#what-is-numerical-differentiation)
  - [Finite Difference Approximations](#finite-difference-approximations)
  - [Error Analysis](#error-analysis)
- [Numerical Differentiation Methods](#-numerical-differentiation-methods)
  - [1. Forward Difference](#1-forward-difference)
  - [2. Backward Difference](#2-backward-difference)
  - [3. Central Difference](#3-central-difference)
- [Method Comparison](#-method-comparison)
- [Applications](#-applications)
- [Implementation Structure](#-implementation-structure)
- [Getting Started](#-getting-started)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

This section of the **Numerical Computing Suite** focuses on **numerical differentiation**—the process of estimating derivatives of functions using discrete data points. When analytical differentiation is difficult or impossible (e.g., for tabulated or experimental data), numerical methods provide practical and efficient solutions.

**The Core Problem**: Given a set of data points `(x₀, y₀), (x₁, y₁), ..., (xₙ, yₙ)`, estimate the derivative $f'(x)$ or higher-order derivatives at specific points.

Numerical differentiation is essential in:
- **Engineering**: Slope, velocity, acceleration from sensor data
- **Physics**: Rate of change in experiments
- **Finance**: Rate of return, volatility
- **Biology**: Growth rates, reaction rates
- **Data Science**: Trend analysis, feature extraction

This collection provides robust finite difference methods for both first and second derivatives, using forward, backward, and central difference formulas.

---

## 🔍 Overview of Numerical Differentiation

Numerical differentiation uses **finite difference formulas** to approximate derivatives based on discrete data. The main approaches are:

| Method              | Data Spacing | Best For                | Accuracy      | Key Feature                |
|---------------------|--------------|------------------------|---------------|----------------------------|
| **Forward**         | Equal        | Near beginning         | $O(h)$        | Uses forward differences   |
| **Backward**        | Equal        | Near end               | $O(h)$        | Uses backward differences  |
| **Central**         | Equal        | Interior points        | $O(h^2)$      | Uses both sides (symmetric)|

---

## 📐 Mathematical Foundation

### What is Numerical Differentiation?

**Numerical differentiation** is the process of estimating the derivative of a function using values at discrete points. When $f(x)$ is known only at certain $x$ values, we use difference quotients to approximate $f'(x)$ and $f''(x)$.

**Formal Definition**: The derivative at $x$ is defined as:

$$
f'(x) = \lim_{h \to 0} \frac{f(x+h) - f(x)}{h}
$$

In practice, with finite $h$, we use difference formulas.

### Finite Difference Approximations

#### **Forward Difference**

$$
f'(x_0) \approx \frac{f(x_1) - f(x_0)}{h}
$$

#### **Backward Difference**

$$
f'(x_n) \approx \frac{f(x_n) - f(x_{n-1})}{h}
$$

#### **Central Difference**

$$
f'(x_i) \approx \frac{f(x_{i+1}) - f(x_{i-1})}{2h}
$$

**Key Concepts**:
- $h$ is the spacing between $x$ values (assumed equal for standard formulas)
- Higher-order differences can be used for better accuracy
- Formulas can be extended for non-uniform spacing (see advanced texts)

### Error Analysis

Numerical differentiation introduces **truncation error** (from Taylor expansion) and **round-off error** (from finite precision arithmetic).

**Error for Forward/Backward Difference**:
$$
f'(x) = \frac{f(x+h) - f(x)}{h} - \frac{h}{2}f''(\xi)
$$
for some $\xi$ in $[x, x+h]$. Error is $O(h)$.

**Error for Central Difference**:
$$
f'(x) = \frac{f(x+h) - f(x-h)}{2h} - \frac{h^2}{6}f'''(\xi)
$$
for some $\xi$ in $[x-h, x+h]$. Error is $O(h^2)$.

**Key Insights**:
- Smaller $h$ reduces truncation error but increases round-off error
- Central difference is generally more accurate
- For noisy data, differentiation amplifies noise—smoothing may be needed

---

## 🛠️ Numerical Differentiation Methods

### 1. Forward Difference


[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-blue?style=for-the-badge)](./First%20and%20Second%20Order%20Derivative%20based%20on%20Forward%20Interpolation/)

#### Theory: Forward Difference Interpolation

Forward difference formulas are derived from the Taylor series expansion about the starting point $x_0$ for equally spaced data ($h = x_{i+1} - x_i$):

$$
f(x_0 + h) = f(x_0) + h f'(x_0) + \frac{h^2}{2!} f''(x_0) + \frac{h^3}{3!} f'''(x_0) + \cdots
$$

By rearranging and solving for derivatives, we obtain:

---

##### 1st Derivative (Forward Difference)

**Formula:**
$$
f'(x_0) \approx \frac{f(x_1) - f(x_0)}{h}
$$

**Higher Accuracy (using more points):**
$$
f'(x_0) \approx \frac{-f(x_2) + 4f(x_1) - 3f(x_0)}{2h}
$$

**Error Analysis:**
- The basic forward difference has a truncation error of $O(h)$:
  $$
  f'(x_0) = \frac{f(x_1) - f(x_0)}{h} - \frac{h}{2} f''(\xi)
  $$
  for some $\xi \in [x_0, x_1]$.
- Using more points (higher-order formula) reduces the error to $O(h^2)$.

##### 2nd Derivative (Forward Difference)

**Formula:**
$$
f''(x_0) \approx \frac{f(x_2) - 2f(x_1) + f(x_0)}{h^2}
$$

**Higher Accuracy (using more points):**
$$
f''(x_0) \approx \frac{-f(x_3) + 4f(x_2) - 5f(x_1) + 2f(x_0)}{h^2}
$$

**Error Analysis:**
- The basic forward difference for the second derivative has a truncation error of $O(h)$:
  $$
  f''(x_0) = \frac{f(x_2) - 2f(x_1) + f(x_0)}{h^2} - \frac{h}{3} f'''(\xi)
  $$
  for some $\xi \in [x_0, x_2]$.
- Higher-order formulas reduce the error to $O(h^2)$ or better.

**Best For:** Estimating derivatives at the beginning of the dataset, especially when only forward data is available.

### 2. Backward Difference


[![View Implementation](https://img.shields.io/badge/📂-View%20Implementation-green?style=for-the-badge)](./First%20and%20Second%20Order%20Derivative%20based%20on%20Backward%20Interpolation/)

#### Theory: Backward Difference Interpolation

Backward difference formulas are derived from the Taylor series expansion about the last point $x_n$ for equally spaced data:

$$
f(x_n - h) = f(x_n) - h f'(x_n) + \frac{h^2}{2!} f''(x_n) - \frac{h^3}{3!} f'''(x_n) + \cdots
$$

By rearranging and solving for derivatives, we obtain:

---

##### 1st Derivative (Backward Difference)

**Formula:**
$$
f'(x_n) \approx \frac{f(x_n) - f(x_{n-1})}{h}
$$

**Higher Accuracy (using more points):**
$$
f'(x_n) \approx \frac{3f(x_n) - 4f(x_{n-1}) + f(x_{n-2})}{2h}
$$

**Error Analysis:**
- The basic backward difference has a truncation error of $O(h)$:
  $$
  f'(x_n) = \frac{f(x_n) - f(x_{n-1})}{h} - \frac{h}{2} f''(\xi)
  $$
  for some $\xi \in [x_{n-1}, x_n]$.
- Using more points (higher-order formula) reduces the error to $O(h^2)$.

##### 2nd Derivative (Backward Difference)

**Formula:**
$$
f''(x_n) \approx \frac{f(x_n) - 2f(x_{n-1}) + f(x_{n-2})}{h^2}
$$

**Higher Accuracy (using more points):**
$$
f''(x_n) \approx \frac{2f(x_n) - 5f(x_{n-1}) + 4f(x_{n-2}) - f(x_{n-3})}{h^2}
$$

**Error Analysis:**
- The basic backward difference for the second derivative has a truncation error of $O(h)$:
  $$
  f''(x_n) = \frac{f(x_n) - 2f(x_{n-1}) + f(x_{n-2})}{h^2} - \frac{h}{3} f'''(\xi)
  $$
  for some $\xi \in [x_{n-2}, x_n]$.
- Higher-order formulas reduce the error to $O(h^2)$ or better.

**Best For:** Estimating derivatives at the end of the dataset, especially when only backward data is available.

### 3. Central Difference

**Formula**:
$$
f'(x_i) \approx \frac{f(x_{i+1}) - f(x_{i-1})}{2h}
$$
$$
f''(x_i) \approx \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2}
$$

**Best For**: Estimating derivatives at interior points. More accurate than forward/backward.

---

## 📊 Method Comparison

| Aspect                | Forward         | Backward        | Central         |
|-----------------------|-----------------|-----------------|-----------------|
| **Spacing**           | Equal           | Equal           | Equal           |
| **Best Location**     | Start           | End             | Middle          |
| **Accuracy**          | $O(h)$          | $O(h)$          | $O(h^2)$        |
| **Noise Sensitivity** | High            | High            | Moderate        |
| **Second Derivative** | Yes             | Yes             | Yes             |
| **Implementation**    | Simple          | Simple          | Simple          |

---

## 🎯 Applications

Numerical differentiation is widely used in:

- **Engineering**: Calculating velocity, acceleration, stress/strain from discrete measurements
- **Physics**: Experimental data analysis, signal processing
- **Finance**: Estimating rates of change in time series
- **Biology**: Growth rates, enzyme kinetics
- **Data Science**: Feature extraction, trend detection
- **Environmental Science**: Rate of change in climate data

---

## 📁 Implementation Structure

Each method folder contains:

```
Method Name/
├── README.md                        # Theory and examples
├── [method]_diff.cpp                # C++ implementation
├── input.txt                        # Test case(s)
├── output.txt                       # Expected output
```

### Common Features

- ✅ File-based I/O for batch processing
- ✅ Multiple test cases per file
- ✅ Output with configurable precision
- ✅ Step-by-step difference table display
- ✅ Error handling for input validation
- ✅ Well-commented, modular code

---

## 🚀 Getting Started

### Prerequisites
- **C++ Compiler**: g++, clang++, or MSVC
- **C++ Standard**: C++11 or later
- **Text Editor/IDE**: VS Code, CLion, or any preferred editor
- **Basic Knowledge**: Calculus and finite differences

### Compilation

Navigate to any method folder and compile:

```bash
# Standard compilation
g++ -o diff [method]_diff.cpp -std=c++11

# With optimization (recommended)
g++ -o diff [method]_diff.cpp -std=c++17 -O2

# With warnings (for development)
g++ -o diff [method]_diff.cpp -std=c++17 -O2 -Wall -Wextra
```

### Execution

```bash
# Run the compiled program
./diff

# On Windows
diff.exe
```

The program:
1. Reads input from `input.txt`
2. Processes all differentiation requests
3. Writes results to `output.txt`
4. Displays results on console

### Input Format

**Basic Format** (`input.txt`):
```
n                   # Number of data points
x₀ y₀              # First data point
x₁ y₁              # Second data point
...                # ...
xₙ₋₁ yₙ₋₁          # Last data point
m                   # Number of points to differentiate
x_diff₁            # First x to differentiate
...                # ...
x_diffₘ            # Last x to differentiate
```

### Output Format

Each run produces:
1. **Header** with method name
2. **Data points table**
3. **Derivative results** for each requested point

### Customization Options

**Adjust Precision**:
```cpp
fout << fixed << setprecision(6);  // Change 6 to desired decimal places
```

**Modify Tolerance**:
```cpp
const double EPSILON = 1e-12;  // Adjust for your numerical requirements
```

---

## 📚 References

- **Numerical Methods for Engineers** by Chapra & Canale
- [Wikipedia: Numerical Differentiation](https://en.wikipedia.org/wiki/Numerical_differentiation)

---

## 👨‍💻 Author

**Abir Hasan Arko**  
Roll: 2207053   
CSE, KUET    
[![GitHub](https://img.shields.io/badge/GitHub-AbirHasanArko-181717?style=flat&logo=github)](https://github.com/AbirHasanArko)

---

<div align="center">

**[⬆ Back to Top](#numerical-differentiation)**

Part of the [Numerical Computing Suite](https://github.com/AbirHasanArko/Numerical-Computing-Suite)

</div>
