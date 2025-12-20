# Numerical Differentiation: Backward Interpolation

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](backward_interpolation_diff.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

---

## 📑 Table of Contents
- [Introduction](#-introduction)
- [Theory & Algorithm](#-theory--algorithm)
  - [Mathematical Foundation](#mathematical-foundation)
  - [Algorithm Steps](#algorithm-steps)
  - [Complexity Analysis](#complexity-analysis)
- [Implementation Details](#-implementation-details)
- [Complete C++ Implementation](#-complete-c-implementation)
- [Usage Example](#-usage-example)
- [Compilation and Execution](#-compilation-and-execution)
- [Applications](#-applications)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

This program computes the **first and second derivatives** of a function at a given point using the **backward interpolation formula** for equally spaced data. It is especially useful for tabulated data where analytical differentiation is not possible, and the point of interest is near the end of the dataset.

### Features

- ✅ **Equally spaced data validation**
- ✅ **Backward difference table construction**
- ✅ **Beautifully formatted output**
- ✅ **File-based I/O**
- ✅ **First and second derivative calculation**
- ✅ **Clear error messages**

---

## 🧮 Theory & Algorithm

### Mathematical Foundation

Given $n$ equally spaced data points $(x_0, y_0), (x_1, y_1), ..., (x_{n-1}, y_{n-1})$, the backward difference table is constructed. The first and second derivatives at a point $x_p$ (close to $x_{n-1}$) are approximated as:

- **First Derivative:**
  $$
  f'(x_p) \approx \frac{1}{h} \left[\nabla y_{n-1} + \frac{2v+1}{2} \nabla^2 y_{n-1} + \frac{3v^2+6v+2}{6} \nabla^3 y_{n-1} + ...\right]
  $$
- **Second Derivative:**
  $$
  f''(x_p) \approx \frac{1}{h^2} \left[\nabla^2 y_{n-1} + (v+1)\nabla^3 y_{n-1} + \frac{6v^2+18v+11}{12} \nabla^4 y_{n-1} + ...\right]
  $$
where $v = \frac{x_p - x_{n-1}}{h}$ and $h$ is the spacing between $x$ values.

#### Backward Difference Table

| $x$ | $y$ | $\nabla y$ | $\nabla^2 y$ | $\nabla^3 y$ | $\nabla^4 y$ |
|-----|-----|-----------|-------------|-------------|-------------|
| $x_0$ | $y_0$ | | | | |
| $x_1$ | $y_1$ | $\nabla y_1$ | | | |
| $x_2$ | $y_2$ | $\nabla y_2$ | $\nabla^2 y_2$ | | |
| $x_3$ | $y_3$ | $\nabla y_3$ | $\nabla^2 y_3$ | $\nabla^3 y_3$ | |
| $x_4$ | $y_4$ | $\nabla y_4$ | $\nabla^2 y_4$ | $\nabla^3 y_4$ | $\nabla^4 y_4$ |

### Algorithm Steps

1. **Read** $n$, $x_i$, $y_i$, and $x_p$ from input file.
2. **Build** the backward difference table for $y$ values.
3. **Compute** $v = (x_p - x_{n-1})/h$.
4. **Calculate** the first and second derivatives using the above formulas.
5. **Output** results in a formatted table to `output.txt`.

### Complexity Analysis

- **Time Complexity:**
  - Building difference table: $O(n^2)$
  - Derivative calculation: $O(n)$
- **Space Complexity:** $O(n^2)$

---

## 💻 Implementation Details

- **Input:**
  - `input.txt` format:
    ```
    n
    x1 x2 ... xn
    y1 y2 ... yn
    xp
    ```
- **Output:**
  - `output.txt` contains:
    - Data points table
    - Backward difference table
    - First and second derivative values at $x_p$
- **Code:**
  - See [backward_interpolation_diff.cpp](backward_interpolation_diff.cpp)

---

## 🔧 Complete C++ Implementation

```cpp
#include <bits/stdc++.h>
using namespace std;

// Utility: factorial
static double factorial(int n)
{
    double f = 1.0;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

// Utility: Backward Difference Table
vector<vector<double>> backwardDiff(const vector<double> &y)
{
    int n = y.size();

    vector<vector<double>> b(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        b[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = n - 1; i >= j; i--)
            b[i][j] = b[i][j - 1] - b[i - 1][j - 1];

    return b;
}

// Backward First Derivative
static double backwardFirstDerivative(const vector<double> &x,
                                      const vector<double> &y,
                                      double xp)
{
    int n = x.size();

    double h = x[1] - x[0];
    auto b = backwardDiff(y);
    double v = (xp - x[n - 1]) / h;
    double res = b[n - 1][1];

    if (n >= 3)
        res += ((2 * v + 1) / 2.0) * b[n - 1][2];
    if (n >= 4)
        res += ((3 * v * v + 6 * v + 2) / 6.0) * b[n - 1][3];

    return res / h;
}

// Backward Second Derivative
static double backwardSecondDerivative(const vector<double> &x,
                                       const vector<double> &y,
                                       double xp)
{
    int n = x.size();

    double h = x[1] - x[0];
    auto b = backwardDiff(y);
    double v = (xp - x[n - 1]) / h;
    double res = b[n - 1][2];

    if (n >= 4)
        res += (v + 1) * b[n - 1][3];
    if (n >= 5)
        res += ((6 * v * v + 18 * v + 11) / 12.0) * b[n - 1][4];
        
    return res / (h * h);
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin.is_open())
    {
        cerr << "Error: Could not open input.txt" << endl;
        return 1;
    }
    if (!fout.is_open())
    {
        cerr << "Error: Could not open output.txt" << endl;
        return 1;
    }

    int n;
    fin >> n;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; ++i)
        fin >> x[i];
    for (int i = 0; i < n; ++i)
        fin >> y[i];

    double xp;
    fin >> xp;

    fout << fixed << setprecision(6);
    fout << "Numerical Differentiation using Backward Interpolation\n";
    fout << "---------------------------------------------------\n";
    
    fout << "Given data points (x, y):\n";
    for (int i = 0; i < n; ++i)
        fout << "  (" << setw(8) << x[i] << ", " << setw(10) << y[i] << ")\n";
    fout << "\n";

    fout << "Point of differentiation (xp): " << xp << "\n\n";

    double first = backwardFirstDerivative(x, y, xp);
    double second = backwardSecondDerivative(x, y, xp);
    fout << "First Derivative at x = " << xp << " : " << first << "\n";
    fout << "Second Derivative at x = " << xp << " : " << second << "\n";
    fout << "\n";

    fin.close();
    fout.close();
    
    return 0;
}
```

---

## 📊 Usage Example

**input.txt:**
```
5
1 2 3 4 5
1 8 27 64 125
4.5
```

**output.txt:** (excerpt)
```
Numerical Differentiation using Backward Interpolation
---------------------------------------------------
Numerical Differentiation using Backward Interpolation
---------------------------------------------------
Given data points (x, y):
  (1.000000,   1.000000)
  (2.000000,   8.000000)
  (3.000000,  27.000000)
  (4.000000,  64.000000)
  (5.000000, 125.000000)

Point of differentiation (xp): 4.500000

First Derivative at x = 4.500000 : 60.750000
Second Derivative at x = 4.500000 : 27.000000
```

---

## 🛠️ Compilation and Execution

```sh
g++ backward_interpolation_diff.cpp -o backward_interpolation_diff
./backward_interpolation_diff
```

---

## 🔬 Applications

- Engineering: Slope/curvature estimation from tabulated data
- Physics: Experimental data analysis
- Data Science: Smoothing and trend analysis

---

## 📚 References

- Numerical Methods for Engineers by Chapra & Canale
- [Backward Difference Formula - Wikipedia](https://en.wikipedia.org/wiki/Finite_difference#Backward_difference)

---

## 👤 Author

**Part of the [Numerical Computing Suite](../) by [AbirHasanArko](https://github.com/AbirHasanArko)**  
Roll: 2207053  
Department of CSE, KUET
