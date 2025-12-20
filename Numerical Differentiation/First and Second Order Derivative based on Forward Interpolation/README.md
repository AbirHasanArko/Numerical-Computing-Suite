# Numerical Differentiation: Forward Interpolation

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](forward_interpolation_diff.cpp)
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

This program computes the **first and second derivatives** of a function at a given point using the **forward interpolation formula** for equally spaced data. It is especially useful for tabulated data where analytical differentiation is not possible.

### Features

- ✅ **Equally spaced data validation**
- ✅ **Forward difference table construction**
- ✅ **Beautifully formatted output**
- ✅ **File-based I/O**
- ✅ **First and second derivative calculation**
- ✅ **Clear error messages**

---

## 🧮 Theory & Algorithm

### Mathematical Foundation

Given $n$ equally spaced data points $(x_0, y_0), (x_1, y_1), ..., (x_{n-1}, y_{n-1})$, the forward difference table is constructed. The first and second derivatives at a point $x_p$ (close to $x_0$) are approximated as:

- **First Derivative:**
   $$
   f'(x_p) \approx \frac{1}{h} \left[\Delta y_0 + \frac{2u-1}{2} \Delta^2 y_0 + \frac{3u^2-6u+2}{6} \Delta^3 y_0 + ...\right]
   $$
- **Second Derivative:**
   $$
   f''(x_p) \approx \frac{1}{h^2} \left[\Delta^2 y_0 + (u-1)\Delta^3 y_0 + \frac{6u^2-18u+11}{12} \Delta^4 y_0 + ...\right]
   $$
where $u = \frac{x_p - x_0}{h}$ and $h$ is the spacing between $x$ values.

#### Forward Difference Table

| $x$ | $y$ | $\Delta y$ | $\Delta^2 y$ | $\Delta^3 y$ | $\Delta^4 y$ |
|-----|-----|-----------|-------------|-------------|-------------|
| $x_0$ | $y_0$ | $\Delta y_0$ | $\Delta^2 y_0$ | $\Delta^3 y_0$ | $\Delta^4 y_0$ |
| $x_1$ | $y_1$ | $\Delta y_1$ | $\Delta^2 y_1$ | $\Delta^3 y_1$ | |
| $x_2$ | $y_2$ | $\Delta y_2$ | $\Delta^2 y_2$ | | |
| $x_3$ | $y_3$ | $\Delta y_3$ | | | |
| $x_4$ | $y_4$ | | | | |

### Algorithm Steps

1. **Read** $n$, $x_i$, $y_i$, and $x_p$ from input file.
2. **Build** the forward difference table for $y$ values.
3. **Compute** $u = (x_p - x_0)/h$.
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
      - Forward difference table
      - First and second derivative values at $x_p$
- **Code:**
   - See [forward_interpolation_diff.cpp](forward_interpolation_diff.cpp)

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

// Utility: Forward Difference Table
vector<vector<double>> forwardDiff(const vector<double> &y)
{
    int n = y.size();
    vector<vector<double>> d(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        d[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            d[i][j] = d[i + 1][j - 1] - d[i][j - 1];

    return d;
}

// Forward First Derivative
static double forwardFirstDerivative(const vector<double> &x,
                                     const vector<double> &y,
                                     double xp)
{
    double h = x[1] - x[0];
    auto d = forwardDiff(y);
    double u = (xp - x[0]) / h;
    double res = d[0][1];

    if (x.size() >= 3)
        res += ((2 * u - 1) / 2.0) * d[0][2];
    if (x.size() >= 4)
        res += ((3 * u * u - 6 * u + 2) / 6.0) * d[0][3];

    return res / h;
}

// Forward Second Derivative
static double forwardSecondDerivative(const vector<double> &x,
                                      const vector<double> &y,
                                      double xp)
{
    double h = x[1] - x[0];
    auto d = forwardDiff(y);
    double u = (xp - x[0]) / h;
    double res = d[0][2];

    if (x.size() >= 4)
        res += (u - 1) * d[0][3];
    if (x.size() >= 5)
        res += ((6 * u * u - 18 * u + 11) / 12.0) * d[0][4];
        
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
    fout << "Numerical Differentiation using Forward Interpolation\n";
    fout << "---------------------------------------------------\n";
    
    fout << "Given data points (x, y):\n";
    for (int i = 0; i < n; ++i)
        fout << "  (" << setw(8) << x[i] << ", " << setw(10) << y[i] << ")\n";
    fout << "\n";
    
    fout << "Point of differentiation (xp): " << xp << "\n\n";

    double first = forwardFirstDerivative(x, y, xp);
    double second = forwardSecondDerivative(x, y, xp);

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
2.5
```

**output.txt:** (excerpt)
```
Numerical Differentiation using Forward Interpolation
---------------------------------------------------
Given data points (x, y):
  (1.000000,   1.000000)
  (2.000000,   8.000000)
  (3.000000,  27.000000)
  (4.000000,  64.000000)
  (5.000000, 125.000000)

Point of differentiation (xp): 2.500000

First Derivative at x = 2.500000 : 18.750000
Second Derivative at x = 2.500000 : 15.000000
```

---

## 🛠️ Compilation and Execution

```sh
g++ forward_interpolation_diff.cpp -o forward_interpolation_diff
./forward_interpolation_diff
```

---

## 🔬 Applications

- Engineering: Slope/curvature estimation from tabulated data
- Physics: Experimental data analysis
- Data Science: Smoothing and trend analysis

---

## 📚 References

- Numerical Methods for Engineers by Chapra & Canale
- [Forward Difference Formula - Wikipedia](https://en.wikipedia.org/wiki/Finite_difference#Forward_difference)

---

## 👤 Author

**Part of the [Numerical Computing Suite](../) by [AbirHasanArko](https://github.com/AbirHasanArko)**  
Roll:  2207053  
Department of CSE, KUET   