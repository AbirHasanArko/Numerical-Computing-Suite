# Secant Method

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](secant-method.cpp)
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

The **Secant Method** is a powerful numerical technique for finding **real roots of nonlinear equations**. It is an **open method** that approximates the derivative using finite differences, making it a derivative-free alternative to Newton-Raphson. The method uses two initial points and iteratively refines the approximation by drawing secants (straight lines) through function values.

This method is particularly effective when:
- The **derivative is difficult or expensive** to compute
- **Faster convergence** than bracketing methods is desired
- Two initial guesses near the root are available
- The function is **smooth and continuous** near the root
- **Superlinear convergence** (faster than linear) is needed

The Secant Method combines the benefits of Newton-Raphson's speed with the simplicity of not requiring derivatives, making it highly practical for engineering and scientific applications.

### Features

- ✅ **Polynomial evaluation** - Efficient computation of polynomial values at any point
- ✅ **Automatic interval detection** - Systematic scanning for sign changes across search range
- ✅ **Derivative-free** - Uses finite difference approximation instead of analytical derivatives
- ✅ **Multiple root finding** - Discovers and computes all real roots within specified range
- ✅ **Superlinear convergence** - Typically faster than False Position and Bisection methods
- ✅ **Cycle detection** - Map-based tracking to detect and handle repeated values
- ✅ **Duplicate root filtering** - Automatic detection and removal of duplicate roots
- ✅ **Division by zero protection** - Safeguards against numerical instabilities
- ✅ **Dual output streams** - Simultaneous output to console and file
- ✅ **High precision formatting** - Configurable decimal precision for scientific accuracy
- ✅ **Robust error handling** - Input validation and meaningful error messages
- ✅ **File-based I/O** - Support for batch processing with structured input files
- ✅ **Clean formatted output** - Professional table formatting with root details

---

## 🧮 Theory & Algorithm

### Mathematical Foundation

The Secant Method is based on **linear interpolation** between two points to approximate the root. Unlike bracketing methods, it does not require the initial points to bracket a root, but it benefits from sign changes for better convergence.

**Core Principle:**

Given two points $x_1$ and $x_2$ with corresponding function values $f(x_1)$ and $f(x_2)$, the method approximates the root by finding where the secant line through these points crosses the x-axis.

**Secant Formula:**

$$x_0 = x_1 - f(x_1) \cdot \frac{x_2 - x_1}{f(x_2) - f(x_1)}$$

This can also be written as:

$$x_0 = \frac{x_1 f(x_2) - x_2 f(x_1)}{f(x_2) - f(x_1)}$$

**Geometric Interpretation:**

The Secant Method draws a straight line (secant) connecting the points $(x_1, f(x_1))$ and $(x_2, f(x_2))$, then finds the x-intercept of this line. This x-intercept becomes the new approximation.

**Derivative Approximation:**

The Secant Method approximates the derivative in Newton-Raphson:

$$f'(x) \approx \frac{f(x_2) - f(x_1)}{x_2 - x_1}$$

Substituting this into Newton's formula $x_{\text{new}} = x_{\text{old}} - \frac{f(x_{\text{old}})}{f'(x_{\text{old}})}$ gives the Secant formula.

**Iteration Process:**

1. Start with two initial points: $x_1$ and $x_2$
2. Calculate new approximation using the Secant formula: $x_0 = x_1 - f(x_1) \cdot \frac{x_2 - x_1}{f(x_2) - f(x_1)}$
3. Check convergence: $|f(x_0)| < \epsilon$
4. Update points: $x_1 = x_2$, $x_2 = x_0$ (shift the two most recent points)
5. Repeat until convergence or cycle detection

**Convergence Criteria:**

The algorithm stops when:

$$|f(x_0)| < \epsilon$$

where:
- $\epsilon = 10^{-6}$ (function value tolerance)

Additionally, the implementation includes **cycle detection** using a map to prevent infinite loops when the method enters a repetitive pattern.

**Comparison with Other Methods:**

| Aspect | Bisection | False Position | Secant | Newton-Raphson |
|--------|-----------|----------------|--------|----------------|
| Derivative Required | No | No | No | Yes |
| Convergence Rate | Linear | Superlinear | Superlinear (≈1.618) | Quadratic (2.0) |
| Initial Points | 2 (bracket) | 2 (bracket) | 2 (any) | 1 |
| Bracketing Required | Yes | Yes | No | No |
| Function Evaluations/Iteration | 1 | 2 | 1 | 1 + derivative |
| Convergence Order | 1.0 | ~1.6 (varies) | φ ≈ 1.618 | 2.0 |
| Best For | Guaranteed slow | Nearly linear | General purpose | When derivative available |

**Convergence Order:**

The Secant Method has a convergence order of $\phi = \frac{1 + \sqrt{5}}{2} \approx 1.618$ (the golden ratio), which means:

$$|e_{n+1}| \approx C |e_n|^{1.618}$$

where $e_n$ is the error at iteration $n$. This is faster than linear (order 1) but slower than quadratic (order 2).

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

3. **Secant Iteration** (for each interval):
   - Initialize: $x_1 = x_i$, $x_2 = x_{i+1}$ (two consecutive points)
   - Create empty visited map for cycle detection
   - While not converged:
     1. Evaluate: $f(x_1) = f(x_1)$, $f(x_2) = f(x_2)$
     2. **Check for division by zero:** If $|f(x_2) - f(x_1)| < 10^{-10}$, break
     3. **Apply Secant formula:** 
        $$x_0 = x_1 - f(x_1) \cdot \frac{x_2 - x_1}{f(x_2) - f(x_1)}$$
     4. Evaluate: $f_0 = f(x_0)$
     5. Round $x_0$ for map: $x_{0,\text{rounded}} = \text{round}(x_0 \times 10^6) / 10^6$
     6. **Convergence/cycle check:**
        - If $|f_0| < \epsilon$ OR $x_{0,\text{rounded}}$ is in visited map:
          - Check for duplicate in roots list
          - Add $x_0$ to roots if unique
          - Break loop
     7. Mark $x_{0,\text{rounded}}$ as visited
     8. **Update for next iteration:** $x_1 = x_2$, $x_2 = x_0$ (shift points)

4. **Root Collection**:
   - Combine roots from direct detection and Secant iterations
   - Filter duplicates (within tolerance $\epsilon$)
   - Sort roots in ascending order

5. **Output Generation**:
   - Display polynomial equation
   - Print final roots table with indices
   - Report completion status

### Complexity Analysis

- **Time Complexity**:
  - Interval scanning: $O(R/h)$ where $R$ is search range, $h$ is step size
  - Single root finding: $O(k)$ where $k$ is number of iterations
  - For $m$ roots: $O(R/h + m \cdot k)$
  - Typically requires **40-60% fewer iterations** than Bisection

- **Space Complexity**: 
  - $O(m + k)$ where $m$ is number of roots
  - $k$ is average size of visited map per root (cycle detection overhead)

- **Convergence Rate**: 
  - **Superlinear** with order $\phi \approx 1.618$
  - Error reduces by factor of approximately $e^{1.618}$ each iteration
  - Significantly faster than Bisection (linear) and False Position
  - Slightly slower than Newton-Raphson (quadratic) but no derivative needed

- **Function Evaluations**:
  - **2 evaluations per iteration** ($f(x_1)$ and $f(x_2)$)
  - More efficient than Newton-Raphson when derivative computation is expensive

- **Accuracy**: 
  - Guaranteed to converge to a root within tolerance $\epsilon$ (when convergent)
  - Final error: $|r_{\text{true}} - r_{\text{computed}}| \lesssim \epsilon$
  - May not converge if initial points are poorly chosen (far from root)

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
       3. Add $c_i \cdot x^{n-i}$ to result using `pow(x, power)`
       4. Return computed value
     - **Complexity:** $O(n)$ where $n$ is polynomial degree

### 2. **Polynomial Display**
   - **`printPolynomial(out, coef)`**
     - Formats and displays polynomial in standard mathematical notation
     - Handles zero coefficients, signs, and special cases
     - **Steps:**
       1. Write "f(x) = "
       2. For each non-zero coefficient:
          - Add sign (+ or -)
          - Add coefficient value (if not ±1 or at constant term)
          - Add variable with power: $x$, $x^2$, $x^3$, etc.
       3. Produce clean output like: $f(x) = x^3 - 6x^2 + 11x - 6$

### 3. **Header Formatting**
   - **`printHeader(out)` / `printHeader()`**
     - Prints styled section headers
     - Overloaded for console and file output
     - Displays "Secant Method" title
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
       4. Handles empty roots list with "No real roots found" message

### 5. **Interval Scanning Engine**
   - **Main scanning loop**
     - Searches entire range for potential root intervals
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
          - Add $x_i$ to intervals list (good starting points)
     - **Output:** List of intervals containing potential roots

### 6. **Secant Core Algorithm**
   - **Main iteration loop for each interval**
     - Implements the Secant Method with cycle detection
     - **Steps:**
       1. Initialize: $x_1 = $ start, $x_2 = $ start + step
       2. Create visited map: `map<f, bool> visited`
       3. While true:
          a. Evaluate: $f(x_1)$ and $f(x_2)$
          b. **Division by zero check:** If $|f(x_2) - f(x_1)| < 10^{-10}$, break
          c. **Apply Secant formula:** 
             $$x_0 = x_1 - f(x_1) \cdot \frac{x_2 - x_1}{f(x_2) - f(x_1)}$$
          d. Evaluate: $f_0 = f(x_0)$
          e. **Round for map comparison:** $x_{0,\text{rounded}} = \text{round}(x_0 \times 10^6) / 10^6$
             - Prevents floating-point precision issues in map lookups
          f. **Convergence/cycle check:**
             - If $|f_0| < \epsilon$: Root found (standard convergence)
             - If $x_{0,\text{rounded}}$ in visited: Cycle detected
             - Add $x_0$ to roots if unique, break loop
          g. Mark $x_{0,\text{rounded}}$ as visited
          h. **Update for next iteration:** 
             - $x_1 = x_2$ (shift: previous $x_2$ becomes new $x_1$)
             - $x_2 = x_0$ (new point becomes new $x_2$)
     - **Key Feature:** Uses two most recent points for next iteration
     - **Cycle Detection Rationale:** 
       - Prevents infinite loops when method oscillates
       - Essential for robust implementation
       - Handles cases where method doesn't converge

### 7. **Cycle Detection Mechanism**
   - **Map-based visited tracking**
     - `map<f, bool> visited`: Stores previously computed values
     - Values rounded to 6 decimal places for comparison
     - **Purpose:**
       - Detect when method enters a repetitive cycle
       - Prevent infinite loops in non-convergent cases
       - Terminate gracefully when convergence stalls
     - **Implementation:**
       ```cpp
       f x0_rounded = round(x0 * 1e6) / 1e6;
       if (visited[x0_rounded]) {
           // Cycle detected, add best approximation and exit
           roots.push_back(x0);
           break;
       }
       visited[x0_rounded] = true;
       ```
     - **Why rounding?** Exact floating-point comparison is unreliable due to precision

### 8. **Duplicate Filtering**
   - **Inline checks throughout code**
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
       - Direct root detection may overlap with Secant results
       - Different starting pairs may find the same root

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
       5. Completion message
     - Ensures professional formatting throughout

### 10. **Data Structures**
   - **`vector<f> coef`**: Stores polynomial coefficients (descending power order)
   - **`vector<f> roots`**: Stores all discovered roots (sorted)
   - **`vector<f> intervals`**: Stores starting points of intervals with sign changes
   - **`map<f, bool> visited`**: Tracks computed values for cycle detection (one per interval)
   - **Type alias:** `#define f double` for flexibility in precision

### 11. **Program Flow**
   1. Display program header to console
   2. Prompt user for input and output filenames
   3. Read input file:
      - Polynomial degree
      - Coefficients
   4. Validate file operations
   5. Write output file header
   6. Write polynomial information to file
   7. Scan search range for:
      - Direct roots (where $f(x) \approx 0$)
      - Intervals with sign changes (good starting points)
   8. For each interval:
      - Apply Secant algorithm with cycle detection
      - Add converged root (if not duplicate)
   9. Sort roots in ascending order
   10. Write roots table to output file
   11. Write completion message
   12. Close files
   13. Display success message to console

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

// Print polynomial in readable form
void printPolynomial(ofstream &out, const vector<f> &coef) {
    int n = coef.size() - 1;
    out << "f(x) = ";

    bool first = true;
    for (int i = 0; i <= n; i++) {
        f c = coef[i];
        int power = n - i;
        if (c == 0) continue;

        if (!first) {
            if (c > 0) out << " + ";
            else out << " - ";
        } else {
            if (c < 0) out << "-";
        }

        if (abs(c) != 1 || power == 0)
            out << abs(c);

        if (power > 0) {
            out << "x";
            if (power > 1)
                out << "^" << power;
        }

        first = false;
    }
    out << "\n\n";
}

// Print header
void printHeader(ofstream &out) {
    out << "============================================\n";
    out << "              Secant Method\n";
    out << "============================================\n\n";
}
void printHeader() {
    cout << "============================================\n";
    cout << "              Secant Method\n";
    cout << "============================================\n\n";
}

// Print roots table
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

int main() {
    string inputFile, outputFile;

    printHeader();
    // File names
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
    for (int i = 0; i <= degree; i++)
        in >> coef[i];

    printHeader(out);

    out << "Polynomial Degree : " << degree << "\n";
    out << "Coefficients      : ";
    for (auto c : coef) out << c << " ";
    out << "\n\n";

    printPolynomial(out, coef);

    // Parameters
    f searchRange = 5000.0, step = 0.5, tolerance = 1e-6;
    vector<f> roots, intervals;

    // Find intervals with sign changes
    for (f i = -searchRange; i <= searchRange; i += step) {
        f fx1 = fun(coef, i);
        f fx2 = fun(coef, i + step);

        // Exact root at grid
        if (abs(fx1) < tolerance) {
            bool duplicate = false;
            for (auto r : roots) if (abs(r - i) < tolerance) { duplicate = true; break; }
            if (!duplicate) roots.push_back(i);
        }
        // Sign change
        else if (fx1 * fx2 < 0.0) {
            intervals.push_back(i);
        }
    }

    // Apply Secant Method
    for (f start : intervals) {
        f x1 = start, x2 = start + step;
        map<f, bool> visited;

        while (true) {
            f fx1 = fun(coef, x1);
            f fx2 = fun(coef, x2);

            // Secant formula
            if (abs(fx2 - fx1) < 1e-10) break; // Avoid division by zero

            f x0 = x1 - fx1 * ((x2 - x1) / (fx2 - fx1));
            f fx0 = fun(coef, x0);

            // Round for map comparison to avoid floating point precision issues
            f x0_rounded = round(x0 * 1e6) / 1e6;

            // Check convergence or repeated value
            if (abs(fx0) < tolerance || visited[x0_rounded]) {
                bool duplicate = false;
                for (auto r : roots) if (abs(r - x0) < tolerance) { duplicate = true; break; }
                if (!duplicate) roots.push_back(x0);
                break;
            }

            visited[x0_rounded] = true;

            // Update for next iteration (Secant method updates)
            x1 = x2;
            x2 = x0;
        }
    }

    sort(roots.begin(), roots.end());

    // Output roots table
    out << "Roots\n";
    out << "============================================\n";
    if (roots.empty()) {
        out << "No real roots found in the given range.\n\n";
    } else {
        printRootsTable(out, roots, coef, tolerance);
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

**Problem:** Find all real roots of the cubic polynomial $f(x) = x^3 - 6x^2 + 11x - 6$ using the Secant Method.

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
./secant
Enter input file name: input1.txt
Enter output file name: output1.txt
```

**Output File (`output1.txt`):**
```
============================================
              Secant Method
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
- ✅ The Secant Method converged rapidly using two-point iterations
- The algorithm identified sign changes in intervals and applied the secant formula iteratively

**Convergence Behavior:**
For this polynomial, the Secant Method typically converges in:
- **Root at x=1:** ~3-4 iterations (superlinear convergence)
- **Root at x=2:** ~3-4 iterations
- **Root at x=3:** ~3-4 iterations

This is significantly faster than Bisection (~20 iterations) and comparable to False Position.

**Verification:**
- $f(1) = 1 - 6 + 11 - 6 = 0$ ✓
- $f(2) = 8 - 24 + 22 - 6 = 0$ ✓
- $f(3) = 27 - 54 + 33 - 6 = 0$ ✓

**Method Efficiency:**
The Secant Method required:
- **No derivative computation** (unlike Newton-Raphson)
- **2 function evaluations per iteration**
- **Superlinear convergence** with order φ ≈ 1.618
- **Automatic cycle detection** prevented any numerical instabilities

---

### Example 2: Polynomial with No Real Roots

**Problem:** Attempt to find real roots of $f(x) = x^2 + 1$ using the Secant Method.

**Mathematical Analysis:**
This polynomial can be factored over complex numbers as:
$$f(x) = (x - i)(x + i)$$

The roots are $x = i$ and $x = -i$ (both complex). There are **no real roots** because:
- $f(x) = x^2 + 1 > 0$ for all real $x$
- The function never crosses the x-axis
- No sign changes exist in any interval

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
./secant
Enter input file name: input2.txt
Enter output file name: output2.txt
```

**Output File (`output2.txt`):**
```
============================================
              Secant Method
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
- The Secant Method cannot be applied without appropriate starting points
- This is mathematically correct since $f(x) = x^2 + 1 \geq 1 > 0$ for all real $x$

**Verification:**
- Discriminant: $\Delta = b^2 - 4ac = 0^2 - 4(1)(1) = -4 < 0$
- Since discriminant is negative, no real roots exist ✓
- Function values: $f(0) = 1$, $f(\pm 1) = 2$, $f(\pm 10) = 101$ (all positive)

**Method Behavior:**
Without sign changes:
- No suitable starting intervals identified
- Secant iteration not initiated
- Clean termination with "No real roots found" message
- Demonstrates robustness of the implementation

---

## 🎯 Compilation and Execution

**Compile:**
```bash
g++ -std=c++17 -O2 secant-method.cpp -o secant
```

**Run:**
```bash
./secant
```

**Alternative (one-liner):**
```bash
g++ -std=c++17 -O2 secant-method.cpp -o secant && ./secant
```

**Compiler Requirements:**
- C++17 or later
- Support for `<bits/stdc++.h>` (commonly available in GCC/MinGW)
- Standard math library (`pow`, `abs`, `fabs`)
- STL containers (`map`, `vector`)

---

## 🔬 Applications

The Secant Method is widely used across various scientific and engineering domains:

1. **Engineering Problems**:
   - Structural analysis (finding critical points without derivatives)
   - Circuit design (solving nonlinear equations)
   - Heat transfer (implicit temperature calculations)
   - Vibration analysis (natural frequency calculations)

2. **Numerical Analysis**:
   - Root-finding when derivatives are unavailable
   - Solving transcendental equations efficiently
   - Inverse function computation
   - Fixed-point problems

3. **Physics Simulations**:
   - Trajectory optimization
   - Energy level calculations
   - Quantum mechanics (Schrödinger equation solutions)
   - Orbital mechanics

4. **Control Systems**:
   - Stability analysis (characteristic equation roots)
   - PID controller tuning
   - Transfer function analysis
   - Root locus calculations

5. **Economics & Finance**:
   - Internal rate of return (IRR) calculations
   - Yield to maturity (bond pricing)
   - Break-even analysis
   - Option pricing models

6. **Chemistry**:
   - Chemical equilibrium (without derivative expressions)
   - Reaction kinetics
   - pH calculations
   - Thermodynamic property estimation

7. **Computer Graphics**:
   - Ray-surface intersection
   - Collision detection
   - Bezier curve analysis
   - Spline evaluation

8. **Aerospace Engineering**:
   - Flight path optimization
   - Aerodynamic coefficient determination
   - Orbit transfer calculations
   - Atmospheric reentry analysis

9. **Optimization Problems**:
   - Finding minima/maxima of complex functions
   - Constraint satisfaction
   - Parameter estimation
   - Model calibration

**Advantages:**
- ✅ **No derivative required** - Only needs function evaluations
- ✅ **Superlinear convergence** - Faster than Bisection and False Position (order ≈ 1.618)
- ✅ **Simple to implement** - Straightforward formula
- ✅ **Efficient** - Fewer iterations than bracketing methods
- ✅ **No bracketing required** - Can start with any two points (though sign change helps)
- ✅ **Practical** - Good balance between speed and simplicity
- ✅ **Robust with cycle detection** - Map-based tracking prevents infinite loops

**Limitations:**
- ❌ **Slower than Newton-Raphson** when derivative is available (order 1.618 vs 2.0)
- ❌ **Cannot find complex roots** - Only detects real roots
- ❌ **May not converge** if initial points are poorly chosen
- ❌ **Division by zero risk** when $f(x_1) \approx f(x_2)$ (mitigated in implementation)
- ❌ **May oscillate or stagnate** for certain functions (handled by cycle detection)
- ❌ **Requires two initial points** - Unlike Newton's one-point start

**When to Use:**
- ✅ Use **Secant Method** when derivative is difficult/expensive to compute
- ✅ Use when you need faster convergence than bracketing methods
- ✅ Use as a derivative-free alternative to Newton-Raphson
- ✅ Use when you have two reasonable starting points near a root
- ❌ Avoid when Newton-Raphson is applicable and derivative is easily computed
- ❌ Avoid for functions with very flat regions (slow convergence)
- ❌ Avoid when guaranteed convergence is critical (use Bisection instead)

**Comparison Summary:**

| Method | Order | Derivative? | Iterations (typical) | Best Use Case |
|--------|-------|-------------|---------------------|---------------|
| **Bisection** | 1.0 | No | ~20 | Guaranteed convergence |
| **False Position** | ~1.5 | No | ~8-10 | Nearly linear functions |
| **Secant** | 1.618 | No | ~4-6 | No derivative available |
| **Newton-Raphson** | 2.0 | Yes | ~3-4 | Derivative available |

The Secant Method offers an excellent **middle ground**: faster than bracketing methods, no derivative needed, and practical for most engineering applications.

---

## 📚 References

- [Secant Method - Wikipedia](https://en.wikipedia.org/wiki/Secant_method)
- Numerical Methods For Engineers by Steven C. Chapra and Raymond P. Canale
- [Root-finding algorithms - Wikipedia](https://en.wikipedia.org/wiki/Root-finding_algorithms)
- Numerical Analysis by Richard L. Burden and J. Douglas Faires
- [Comparison of root-finding methods](https://en.wikipedia.org/wiki/Root-finding_algorithm#Comparison)

---

## 👤 Author

**Part of the [Numerical Computing Suite](../) by [MD. Shomik Shahriar](https://github.com/Hapi-Guy)**  
Roll: 2207041  
Department of CSE, KUET
