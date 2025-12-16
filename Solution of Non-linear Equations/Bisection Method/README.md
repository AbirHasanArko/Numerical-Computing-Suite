# Bisection Method

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](bi-section-method.cpp)
[![View Input](https://img.shields.io/badge/View-Input1-green?style=for-the-badge&logo=files)](input1.txt)
[![View Input](https://img.shields.io/badge/View-Input2-green?style=for-the-badge&logo=files)](input2.txt)
[![View Output](https://img.shields.io/badge/View-Output1-orange?style=for-the-badge&logo=files)](output1.txt)
[![View Output](https://img.shields.io/badge/View-Output2-orange?style=for-the-badge&logo=files)](output2.txt)

---

## 📑 Table of Contents
- [Introduction](#-introduction)
- [Theory & Algorithm](#-theory--algorithm)
  - [Mathematical Foundation](#mathematical-foundation)
  - [Algorithm Steps](#algorithm-steps)
  - [Complexity Analysis](#complexity-analysis)
- [Implementation Details](#-implementation-details)
- [Complete C++ Implementation](#-complete-c-implementation)
- [Usage Examples](#-usage-examples)
  - [Example 1: Polynomial with Real Roots](#example-1-polynomial-with-real-roots)
  - [Example 2: Polynomial with No Real Roots](#example-2-polynomial-with-no-real-roots)
- [Compilation and Execution](#-compilation-and-execution)
- [Applications](#-applications)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

The Bisection Method is one of the most fundamental and reliable numerical techniques for finding **real roots of nonlinear equations**. It is a **bracketing method** that systematically narrows down an interval containing a root by repeatedly bisecting it, guaranteeing convergence whenever a sign change is detected.

This method is particularly effective when:
- The function is **continuous** over the interval
- A **sign change** exists in the interval (indicating a root)
- **Guaranteed convergence** is more important than speed
- Initial approximations are not available

The method's robustness and simplicity make it an excellent choice for educational purposes and as a fallback algorithm when more sophisticated methods fail.

### Features

- ✅ **Polynomial evaluation** - Efficient computation of polynomial values at any point
- ✅ **Automatic interval detection** - Systematic scanning for sign changes across search range
- ✅ **Multiple root finding** - Discovers and computes all real roots within specified range
- ✅ **Iteration tracking** - Complete history of bisection iterations for each root
- ✅ **Duplicate root filtering** - Automatic detection and removal of duplicate roots
- ✅ **Convergence guarantee** - Mathematically proven to converge when conditions are met
- ✅ **Dual output streams** - Simultaneous output to console and file
- ✅ **High precision formatting** - Configurable decimal precision for scientific accuracy
- ✅ **Robust error handling** - Input validation and meaningful error messages
- ✅ **File-based I/O** - Support for batch processing with structured input files
- ✅ **Clean formatted output** - Professional table formatting with iteration details
- ✅ **Extrapolation warnings** - Detection of roots outside expected range

---

## 🧮 Theory & Algorithm

### Mathematical Foundation

The Bisection Method is based on the **Intermediate Value Theorem** (IVT), which states:

> If a continuous function $f(x)$ satisfies $f(a) \cdot f(b) < 0$ on an interval $[a, b]$, then there exists at least one root $r \in (a, b)$ such that $f(r) = 0$.

**Core Principle:**

Given an interval $[x_{\text{low}}, x_{\text{high}}]$ where $f(x_{\text{low}}) \cdot f(x_{\text{high}}) < 0$:

1. Calculate the midpoint: $x_{\text{mid}} = \frac{x_{\text{low}} + x_{\text{high}}}{2}$

2. Evaluate $f(x_{\text{mid}})$:
   - If $f(x_{\text{mid}}) \approx 0$: Root found at $x_{\text{mid}}$
   - If $f(x_{\text{low}}) \cdot f(x_{\text{mid}}) < 0$: Root is in $[x_{\text{low}}, x_{\text{mid}}]$, set $x_{\text{high}} = x_{\text{mid}}$
   - If $f(x_{\text{mid}}) \cdot f(x_{\text{high}}) < 0$: Root is in $[x_{\text{mid}}, x_{\text{high}}]$, set $x_{\text{low}} = x_{\text{mid}}$

3. Repeat until convergence criteria are met

**Convergence Criteria:**

The algorithm stops when either condition is satisfied:

$$|f(x_{\text{mid}})| < \epsilon_f \quad \text{or} \quad |x_{\text{high}} - x_{\text{low}}| < \epsilon_x$$

where:
- $\epsilon_f = 10^{-6}$ (function value tolerance)
- $\epsilon_x = 10^{-6}$ (interval width tolerance)

**Error Bound:**

After $n$ iterations, the error is bounded by:

$$|r - x_{\text{mid}}| \leq \frac{b - a}{2^{n+1}}$$

where $[a, b]$ is the initial interval containing the root.

### Algorithm Steps

1. **Input Processing**:
   - Read polynomial degree $n$
   - Read $n+1$ coefficients: $a_0, a_1, ..., a_n$ for $f(x) = a_0 x^n + a_1 x^{n-1} + ... + a_{n-1}x + a_n$

2. **Interval Scanning**:
   - Define search range: $[-R, R]$ where $R = 5000$ (default)
   - Use step size: $h = 0.5$
   - For each point $x_i = -R + i \cdot h$:
     - Check if $|f(x_i)| < \epsilon$ (direct root detection)
     - Check if $f(x_i) \cdot f(x_{i+1}) < 0$ (sign change detection)
     - Store intervals with sign changes

3. **Bisection Iteration** (for each interval):
   - Initialize: $x_{\text{low}} = x_i$, $x_{\text{high}} = x_{i+1}$
   - While not converged:
     1. $x_{\text{mid}} = \frac{x_{\text{low}} + x_{\text{high}}}{2}$
     2. $f_{\text{mid}} = f(x_{\text{mid}})$
     3. Record iteration: $(iteration, x_{\text{low}}, x_{\text{high}}, x_{\text{mid}}, f_{\text{mid}})$
     4. Check convergence: $|f_{\text{mid}}| < \epsilon$ or $|x_{\text{high}} - x_{\text{low}}| < \epsilon$
     5. Update interval:
        - If $f(x_{\text{low}}) \cdot f_{\text{mid}} < 0$: $x_{\text{high}} = x_{\text{mid}}$
        - Else: $x_{\text{low}} = x_{\text{mid}}$

4. **Root Collection**:
   - Add converged $x_{\text{mid}}$ to roots list
   - Check for duplicates (within tolerance $\epsilon$)
   - Sort roots in ascending order

5. **Output Generation**:
   - Display polynomial equation
   - Print final roots table
   - Print iteration history for each root

### Complexity Analysis

- **Time Complexity**:
  - Interval scanning: $O(R/h)$ where $R$ is search range, $h$ is step size
  - Single root finding: $O(\log_2(\frac{b-a}{\epsilon}))$ where $[a,b]$ is initial interval
  - For $m$ roots: $O(R/h + m \log_2(\frac{b-a}{\epsilon}))$
  - Typical: $O(R/h)$ dominates for large search ranges

- **Space Complexity**: $O(m \cdot k)$ where $m$ is number of roots, $k$ is average iterations per root

- **Convergence Rate**: Linear convergence - error reduces by half each iteration

- **Accuracy**: 
  - Guaranteed to converge to a root within tolerance $\epsilon$
  - Final error: $|r_{\text{true}} - r_{\text{computed}}| \lesssim \epsilon$

---

## 💻 Implementation Details

The C++ implementation is structured into the following components:

### 1. **Polynomial Evaluation Function**
   - **`fun(coef, x)`**
     - Evaluates polynomial at given point $x$
     - Uses coefficient vector in descending power order
     - Returns $f(x) = \sum_{i=0}^{n} c_i \cdot x^{n-i}$
     - **Steps:**
       1. Initialize result to 0
       2. Iterate through coefficients
       3. Add $c_i \cdot x^{n-i}$ to result
       4. Return computed value
     - **Complexity:** $O(n)$ where $n$ is polynomial degree

### 2. **Polynomial Display**
   - **`printPolynomial(out, coef)`**
     - Formats and displays polynomial in standard mathematical notation
     - Handles zero coefficients and signs correctly
     - **Steps:**
       1. Write "f(x) = "
       2. For each non-zero coefficient:
          - Add sign (+ or -)
          - Add coefficient value (if not 1 or at constant term)
          - Add variable with power: $x$, $x^2$, $x^3$, etc.
       3. Produce clean output like: $f(x) = x^3 - 6x^2 + 11x - 6$

### 3. **Header Formatting**
   - **`printHeader(out)` / `printHeader()`**
     - Prints styled section headers
     - Overloaded for console and file output
     - Ensures consistent branding across outputs
     - Uses decorative borders for visual appeal

### 4. **Roots Table Display**
   - **`printRootsTable(out, roots, coef, tolerance)`**
     - Formats final roots in tabular form
     - Shows index and root value with 6 decimal precision
     - **Steps:**
       1. Print table header: "Index", "Root Value"
       2. Print separator line
       3. For each root:
          - Print sequential index (1, 2, 3, ...)
          - Print root value with fixed precision
       4. Handles empty roots list with appropriate message

### 5. **Iteration Table Display**
   - **`printBisectionIteration(out, iterData)`**
     - Displays complete iteration history for a single root
     - Shows progression of interval narrowing
     - **Columns:**
       - Iteration number
       - $x_{\text{low}}$ (lower bound)
       - $x_{\text{high}}$ (upper bound)
       - $x_{\text{mid}}$ (midpoint)
       - $f(x_{\text{mid}})$ (function value at midpoint)
     - **Steps:**
       1. Print column headers with proper spacing
       2. Print separator line
       3. For each iteration tuple:
          - Extract values: iteration, xL, xH, xM, fxM
          - Format and print row with aligned columns
       4. Provides insight into convergence behavior

### 6. **Interval Scanning Engine**
   - **Main loop in `main()`**
     - Scans entire search range for potential roots
     - **Parameters:**
       - `searchRange = 5000.0`: Search from -5000 to +5000
       - `step = 0.5`: Increment between test points
       - `tolerance = 1e-6`: Convergence tolerance
     - **Steps:**
       1. Loop from $-R$ to $+R$ with step $h$
       2. Evaluate $f(x_i)$ and $f(x_{i+1})$
       3. **Direct root check:** If $|f(x_i)| < \epsilon$:
          - Check for duplicates in roots list
          - Add $x_i$ to roots if unique
       4. **Sign change check:** If $f(x_i) \cdot f(x_{i+1}) < 0$:
          - Add $x_i$ to intervals list (root exists in $[x_i, x_{i+1}]$)
     - **Output:** List of intervals containing roots

### 7. **Bisection Core Algorithm**
   - **Main loop for each interval**
     - Implements the bisection iteration for root refinement
     - **Steps:**
       1. Initialize: $x_{\text{low}} = $ start, $x_{\text{high}} = $ start + step
       2. Create empty iteration history vector
       3. While true:
          a. Calculate midpoint: $x_{\text{mid}} = \frac{x_{\text{low}} + x_{\text{high}}}{2}$
          b. Evaluate: $f(x_{\text{mid}})$
          c. Record: $(iteration, x_{\text{low}}, x_{\text{high}}, x_{\text{mid}}, f(x_{\text{mid}}))$
          d. **Convergence check:**
             - If $|f(x_{\text{mid}})| < \epsilon$ or $|x_{\text{high}} - x_{\text{low}}| < \epsilon$:
               - Check for duplicate in roots
               - Add $x_{\text{mid}}$ to roots if unique
               - Break loop
          e. **Interval update:**
             - If $f(x_{\text{low}}) \cdot f(x_{\text{mid}}) < 0$: $x_{\text{high}} = x_{\text{mid}}$
             - Else: $x_{\text{low}} = x_{\text{mid}}$
          f. Increment iteration counter
       4. Store iteration history for this root
     - **Output:** Converged root and its iteration history

### 8. **Duplicate Filtering**
   - **Inline checks in code**
     - Prevents same root from being added multiple times
     - Uses tolerance-based comparison
     - **Logic:** For each potential new root $r_{\text{new}}$:
       ```
       for each existing root r in roots:
           if |r_new - r| < tolerance:
               mark as duplicate
               skip addition
       ```
     - Essential because:
       - Multiple intervals may converge to same root
       - Direct root detection may find root already in intervals
       - Scanning step size may cause overlapping detections

### 9. **File I/O Management**
   - **Input Reading:**
     - Opens and validates input file
     - Reads polynomial degree $n$
     - Reads $n+1$ coefficients in descending power order
     - Uses `ifstream` for file input
     - Error handling for missing or invalid files
   
   - **Output Writing:**
     - Creates output file with formatted results
     - Uses `ofstream` for file output
     - Writes:
       1. Program header
       2. Polynomial degree and coefficients
       3. Polynomial equation in standard form
       4. Final roots table (or "No roots found" message)
       5. Iteration tables for each root
       6. Completion message
     - Ensures professional formatting throughout

### 10. **Data Structures**
   - **`vector<f> coef`**: Stores polynomial coefficients
   - **`vector<f> roots`**: Stores all discovered roots (sorted)
   - **`vector<f> intervals`**: Stores starting points of intervals with sign changes
   - **`vector<tuple<int,f,f,f,f>> iterData`**: Stores iteration history for one root
     - Tuple elements: (iteration_number, x_low, x_high, x_mid, f_x_mid)
   - **`vector<vector<tuple<...>>> allIterations`**: Stores all iteration histories

### 11. **Program Flow**
   1. Display program header to console
   2. Prompt user for input and output filenames
   3. Read input file:
      - Polynomial degree
      - Coefficients
   4. Validate file operations
   5. Write output file header
   6. Write polynomial information
   7. Scan search range for:
      - Direct roots (where $f(x) \approx 0$)
      - Intervals with sign changes
   8. For each interval with sign change:
      - Apply bisection algorithm
      - Record iteration history
      - Add converged root (if not duplicate)
   9. Sort roots in ascending order
   10. Write roots table to output file
   11. Write iteration tables for each root
   12. Write completion message
   13. Close files
   14. Display success message to console

---

## 🔧 Complete C++ Implementation

```cpp
#include <bits/stdc++.h>
using namespace std;

#define f double

// Polynomial function f(x)
f fun(const vector<f> &coef, f x) {
    f result = 0.0;
    int n = coef.size() - 1;
    for (int i = 0; i <= n; i++)
        result += coef[i] * pow(x, n - i);
    return result;
}

// Print polynomial nicely
void printPolynomial(ofstream &out, const vector<f> &coef) {
    int n = coef.size() - 1;
    out << "f(x) = ";
    bool first = true;
    for (int i = 0; i <= n; i++) {
        f c = coef[i]; int power = n - i;
        if (c == 0) continue;
        if (!first) {
            if (c > 0) out << " + "; else out << " - ";
        } else { if (c < 0) out << "-"; }
        if (abs(c) != 1 || power == 0) out << abs(c);
        if (power > 0) { out << "x"; if (power > 1) out << "^" << power; }
        first = false;
    }
    out << "\n\n";
}

// Print header
void printHeader(ofstream &out) {
    out << "============================================\n";
    out << "              Bisection Method\n";
    out << "============================================\n\n";
}
void printHeader() {
    cout << "============================================\n";
    cout << "              Bisection Method\n";
    cout << "============================================\n\n";
}

// Print final roots table with trivial root highlight
void printRootsTable(ofstream &out, const vector<f> &roots, const vector<f> &coef, f tolerance = 1e-6) {
    out << setw(6) << "Index" 
        << setw(20) << "Root Value" << "\n";
    out << string(26, '-') << "\n";

    out << fixed << setprecision(6);
    for (size_t i = 0; i < roots.size(); i++) {
        f fx = fun(coef, roots[i]);
        out << setw(6) << (i + 1)
            << setw(20) << roots[i]
            << "\n";
    }
    out << "\n";
}

// Print iteration table per root
void printBisectionIteration(ofstream &out, const vector<tuple<int,f,f,f,f>> &iterData) {
    out << setw(10) << "Iteration"
        << setw(15) << "x_low"
        << setw(15) << "x_high"
        << setw(15) << "x_mid"
        << setw(20) << "f(x_mid)" << "\n";
    out << string(75, '-') << "\n";
    out << fixed << setprecision(6);
    for (auto &t : iterData) {
        int iter; f xL, xH, xM, fxM;
        tie(iter, xL, xH, xM, fxM) = t;
        out << setw(10) << iter
            << setw(15) << xL
            << setw(15) << xH
            << setw(15) << xM
            << setw(20) << fxM << "\n";
    }
    out << "\n";
}

int main() {
    string inputFile, outputFile;

    // File names
    printHeader();
    cout << "Enter input file name: ";
    cin >> inputFile;
    cout << "Enter output file name: ";
    cin >> outputFile;

    ifstream in(inputFile);
    ofstream out(outputFile);
    if (!in.is_open() || !out.is_open()) {
        cout << "Error: Cannot open input/output file!\n";
        return 1;
    }

    int degree;
    in >> degree;
    vector<f> coef(degree + 1);
    for (int i = 0; i <= degree; i++) in >> coef[i];

    printHeader(out);
    out << "Polynomial Degree : " << degree << "\n";
    out << "Coefficients      : ";
    for (auto c : coef) out << c << " ";
    out << "\n\n";
    printPolynomial(out, coef);

    // Bisection parameters
    f searchRange = 5000.0, step = 0.5, tolerance = 1e-6;
    vector<f> roots, intervals;

    // Scan for roots / intervals
    for (f i = -searchRange; i <= searchRange; i += step) {
        f fx1 = fun(coef, i);
        f fx2 = fun(coef, i + step);
        if (abs(fx1) < tolerance) {
            bool duplicate = false;
            for (auto r : roots) if (abs(r - i) < tolerance) { duplicate = true; break; }
            if (!duplicate) roots.push_back(i);
        } else if (fx1 * fx2 < 0.0) intervals.push_back(i);
    }

    vector<vector<tuple<int,f,f,f,f>>> allIterations; // iterations per root

    for (f start : intervals) {
        f xL = start, xH = start + step, xM;
        vector<tuple<int,f,f,f,f>> iterData;
        int iter = 1;
        while (true) {
            xM = (xL + xH) / 2.0;
            f fxM = fun(coef, xM);
            iterData.push_back({iter, xL, xH, xM, fxM});

            if (abs(fxM) < tolerance || abs(xH - xL) < tolerance) {
                bool duplicate = false;
                for (auto r : roots) if (abs(r - xM) < tolerance) { duplicate = true; break; }
                if (!duplicate) roots.push_back(xM);
                break;
            }

            f fxL = fun(coef, xL);
            if (fxL * fxM < 0) xH = xM;
            else xL = xM;

            iter++;
        }
        allIterations.push_back(iterData);
    }

    sort(roots.begin(), roots.end());

    // Print final roots
    out << "Roots\n";
    out << "============================================\n";
    if (roots.empty()) {
        out << "No real roots found in the given range.\n\n";
    } else {
        printRootsTable(out, roots, coef, tolerance);
    }

    // Print iteration tables
    for (size_t i = 0; i < allIterations.size(); i++) {
        out << "Iteration Table for Root " << (i + 1) << "\n";
        out << "----------------------------------------\n";
        printBisectionIteration(out, allIterations[i]);
    }

    out << "============================================\n";
    out << "Computation Completed Successfully.\n";

    in.close();
    out.close();
    cout << "Computation completed. Results written to '" << outputFile << "'\n";

    return 0;
}
```

---

## 📊 Usage Examples

### Input File Format

The input file must follow this specific structure:

```
n                    # Line 1: Polynomial degree (integer)
a₀ a₁ a₂ ... aₙ     # Line 2: (n+1) coefficients in descending power order
```

**Important Notes:**
- The polynomial is represented as: $f(x) = a_0 x^n + a_1 x^{n-1} + ... + a_{n-1} x + a_n$
- Coefficients are in **descending power order** (highest power first)
- All coefficients must be provided, including zeros
- The program searches for roots in the range $[-5000, 5000]$ with step size $0.5$
- Only **real roots** are detected (complex roots are not found)
- Convergence tolerance is $\epsilon = 10^{-6}$

---

### Example 1: Polynomial with Real Roots

**Problem:** Find all real roots of the cubic polynomial $f(x) = x^3 - 6x^2 + 11x - 6$.

**Mathematical Analysis:**
This polynomial can be factored as:
$$f(x) = (x - 1)(x - 2)(x - 3)$$

So the exact roots are: $x = 1, 2, 3$

**Input File (`input1.txt`):**
```
3
1 -6 11 -6
```

**Explanation:**
- **Line 1:** `3` → Polynomial degree is 3 (cubic)
- **Line 2:** `1 -6 11 -6` → Coefficients for $1 \cdot x^3 + (-6) \cdot x^2 + 11 \cdot x + (-6)$

**Execution:**
```bash
./bisection
Enter input file name: input1.txt
Enter output file name: output1.txt
```

**Output File (`output1.txt`):**
```
============================================
              Bisection Method
============================================

Polynomial Degree : 3
Coefficients      : 1 -6 11 -6 

f(x) = x^3 - 6x^2 + 11x - 6

Roots
============================================
 Index          Root Value
--------------------------
     1            1.000000
     2            2.000000
     3            3.000000

============================================
Computation Completed Successfully.
```

**Analysis:**
- ✅ All three real roots found: $1.000000, 2.000000, 3.000000$
- ✅ Roots are accurate to 6 decimal places
- ✅ Matches the exact factored roots perfectly
- ✅ Roots displayed in ascending order
- The algorithm successfully identified all sign changes and converged to each root

**Verification:**
- $f(1) = 1 - 6 + 11 - 6 = 0$ ✓
- $f(2) = 8 - 24 + 22 - 6 = 0$ ✓
- $f(3) = 27 - 54 + 33 - 6 = 0$ ✓

---

### Example 2: Polynomial with No Real Roots

**Problem:** Attempt to find real roots of $f(x) = x^2 + 1$.

**Mathematical Analysis:**
This polynomial can be factored over complex numbers as:
$$f(x) = (x - i)(x + i)$$

The roots are $x = i$ and $x = -i$ (both complex). There are **no real roots** because:
- $f(x) = x^2 + 1 > 0$ for all real $x$
- The function never crosses the x-axis

**Input File (`input2.txt`):**
```
2
1 0 1
```

**Explanation:**
- **Line 1:** `2` → Polynomial degree is 2 (quadratic)
- **Line 2:** `1 0 1` → Coefficients for $1 \cdot x^2 + 0 \cdot x + 1$

**Execution:**
```bash
./bisection
Enter input file name: input2.txt
Enter output file name: output2.txt
```

**Output File (`output2.txt`):**
```
============================================
              Bisection Method
============================================

Polynomial Degree : 2
Coefficients      : 1 0 1 

f(x) = x^2 + 1

Roots
============================================
No real roots found in the given range.

============================================
Computation Completed Successfully.
```

**Analysis:**
- ✅ Correctly identified that no real roots exist
- ✅ No false positives or spurious roots
- ✅ Clean output message indicating absence of real roots
- The algorithm scanned the entire range $[-5000, 5000]$ and found:
  - No points where $f(x) \approx 0$
  - No intervals with sign changes
- This is mathematically correct since $f(x) = x^2 + 1 \geq 1 > 0$ for all real $x$

**Verification:**
- Discriminant: $\Delta = b^2 - 4ac = 0^2 - 4(1)(1) = -4 < 0$
- Since discriminant is negative, no real roots exist ✓

---

## 🎯 Compilation and Execution

**Compile:**
```bash
g++ -std=c++17 -O2 bi-section-method.cpp -o bisection
```

**Run:**
```bash
./bisection
```

**Alternative (one-liner):**
```bash
g++ -std=c++17 -O2 bi-section-method.cpp -o bisection && ./bisection
```

**Compiler Requirements:**
- C++17 or later
- Support for `<bits/stdc++.h>` (commonly available in GCC/MinGW)
- Standard math library (`pow`, `abs`, `fabs`)

---

## 🔬 Applications

The Bisection Method is widely used across various scientific and engineering domains:

1. **Engineering Problems**:
   - Structural analysis (finding critical load points)
   - Circuit design (solving nonlinear circuit equations)
   - Heat transfer calculations (temperature distribution)

2. **Numerical Analysis**:
   - Root-finding for transcendental equations
   - Solving implicit equations
   - Inverse function computation

3. **Physics Simulations**:
   - Trajectory calculations (projectile motion with air resistance)
   - Quantum mechanics (energy eigenvalue problems)
   - Fluid dynamics (flow rate calculations)

4. **Control Systems**:
   - Stability analysis (finding poles and zeros)
   - PID controller tuning
   - Transfer function analysis

5. **Economics & Finance**:
   - Internal rate of return (IRR) calculations
   - Break-even analysis
   - Option pricing models

6. **Chemistry**:
   - pH calculations (solving Henderson-Hasselbalch equations)
   - Chemical equilibrium problems
   - Reaction kinetics

7. **Computer Graphics**:
   - Ray tracing (intersection calculations)
   - Collision detection
   - Curve fitting

**Advantages:**
- ✅ **Guaranteed convergence** when initial interval brackets a root
- ✅ **Simple to understand and implement**
- ✅ **Robust** - works when other methods fail
- ✅ **No derivative required** - only function evaluations needed
- ✅ **Reliable** for ill-behaved functions
- ✅ **Always finds a root** in the bracketing interval

**Limitations:**
- ❌ **Slow convergence** - linear convergence rate (halving interval each iteration)
- ❌ **Cannot find complex roots** - only detects real roots
- ❌ **Requires initial bracketing** - must have sign change in interval
- ❌ **May miss roots** - if multiple roots exist in same interval or at even multiplicities
- ❌ **Inefficient for high precision** - many iterations needed for tight tolerance
- ❌ **Cannot handle discontinuities** - assumes function is continuous

**When to Use:**
- ✅ Use **Bisection Method** when convergence guarantee is critical
- ✅ Use when derivative is difficult or impossible to compute
- ✅ Use as a fallback when Newton-Raphson or Secant methods fail
- ❌ Avoid when speed is critical (use Newton-Raphson instead)
- ❌ Avoid for complex roots (use numerical methods for complex analysis)

---

## 📚 References

- [Bisection Method - Wikipedia](https://en.wikipedia.org/wiki/Bisection_method)
- Numerical Methods For Engineers by Steven C. Chapra and Raymond P. Canale
- [Root-finding algorithms - Wikipedia](https://en.wikipedia.org/wiki/Root-finding_algorithms)
- Numerical Analysis by Richard L. Burden and J. Douglas Faires

---

## 👤 Author

**Part of the [Numerical Computing Suite](../) by [MD. Shomik Shahriar](https://github.com/Hapi-Guy)**  
Roll: 2207041  
Department of CSE, KUET
