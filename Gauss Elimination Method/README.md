# Gauss Elimination Method

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](gauss-elimination-method.cpp)
[![View Input](https://img.shields.io/badge/View-Input-green?style=for-the-badge&logo=files)](input.txt)
[![View Output](https://img.shields.io/badge/View-Output-orange?style=for-the-badge&logo=files)](output.txt)

## 📖 Introduction

The **Gauss Elimination Method** (also known as **Gaussian Elimination**) is a fundamental algorithm in linear algebra used to solve systems of linear equations. This method transforms the coefficient matrix into an upper triangular form through a series of row operations, followed by back substitution to find the solution.

This implementation in C++ solves systems of linear equations and can detect three types of solutions:
- **Unique Solution** - The system has exactly one solution
- **No Solution** - The system is inconsistent
- **Infinite Solutions** - The system has infinitely many solutions

## 🧮 Theory & Algorithm

### Mathematical Foundation

Given a system of linear equations:
```
a₁₁x₁ + a₁₂x₂ + ... + a₁ₙxₙ = b₁
a₂₁x₁ + a₂₂x₂ + ... + a₂ₙxₙ = b₂
... 
aₙ₁x₁ + aₙ₂x₂ + ... + aₙₙxₙ = bₙ
```

The method converts the augmented matrix [A|b] into row echelon form using elementary row operations.

### Algorithm Steps

1. **Forward Elimination Phase:**
   - For each pivot position (i, i):
     - **Partial Pivoting**: Find the row with the maximum absolute value in column i (from row i to n)
     - Swap the current row with the maximum row
     - Eliminate all entries below the pivot by subtracting multiples of the pivot row
     - Calculate the multiplication factor: `factor = a[k][i] / a[i][i]`
     - Update row k: `a[k][j] = a[k][j] - factor × a[i][j]`

2. **Solution Detection Phase:**
   - Calculate the rank of the coefficient matrix
   - Check for inconsistency: If any row has all zeros in coefficients but non-zero constant → **No Solution**
   - Check for infinite solutions: If rank < n → **Infinite Solutions**
   - Otherwise → **Unique Solution**, proceed to back substitution

3. **Back Substitution Phase** (for unique solutions):
   - Start from the last equation: `xₙ = aₙₙ / bₙ`
   - For each variable from n-1 to 1:
     - `xᵢ = (bᵢ - Σ(aᵢⱼ × xⱼ)) / aᵢᵢ` where j goes from i+1 to n

### Time Complexity
- **Forward Elimination**: O(n³)
- **Back Substitution**: O(n²)
- **Overall**: O(n³)

### Space Complexity
- O(n²) for storing the augmented matrix

## 💻 Code Explanation

### Key Components

```cpp
vector<vector<double>> a(n, vector<double>(n + 1));
```
- Creates an augmented matrix of size n×(n+1) to store coefficients and constants

### Partial Pivoting (Lines 44-48)
```cpp
int maxRow = i;
for (int k = i + 1; k < n; k++)
    if (fabs(a[k][i]) > fabs(a[maxRow][i]))
        maxRow = k;
swap(a[i], a[maxRow]);
```
- Finds the row with the largest absolute value in the current column
- Swaps it with the current pivot row
- This improves numerical stability and avoids division by very small numbers

### Forward Elimination (Lines 53-58)
```cpp
for (int k = i + 1; k < n; k++)
{
    double factor = a[k][i] / a[i][i];
    for (int j = i; j <= n; j++)
        a[k][j] -= factor * a[i][j];
}
```
- Eliminates elements below the pivot
- Creates zeros in column i for all rows below row i

### Solution Detection (Lines 74-93)
- Counts the rank by checking non-zero rows
- Detects inconsistency: `0 = non-zero constant`
- Detects infinite solutions when rank < n

### Back Substitution (Lines 103-109)
```cpp
for (int i = n - 1; i >= 0; i--)
{
    x[i] = a[i][n];
    for (int j = i + 1; j < n; j++)
        x[i] -= a[i][j] * x[j];
    x[i] /= a[i][i];
}
```
- Solves for variables from bottom to top
- Substitutes already-known variables into equations above

## 🧪 Test Cases

### Test Case 1: Unique Solution
**Input:**
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

**System of Equations:**
```
2x₁ + x₂ - x₃ = 8
-3x₁ - x₂ + 2x₃ = -11
-2x₁ + x₂ + 2x₃ = -3
```

**Output:**
```
x₁ = 2.00
x₂ = 3.00
x₃ = -1. 00
```

### Test Case 2: No Solution (Inconsistent System)
**Input:**
```
3
1 1 1 2
2 2 2 4
1 1 1 5
```

**System of Equations:**
```
x₁ + x₂ + x₃ = 2
2x₁ + 2x₂ + 2x₃ = 4
x₁ + x₂ + x₃ = 5
```

**Output:**
```
No Solution
```

**Explanation:** The first and third equations contradict each other (both have identical coefficients but different constants). 

### Test Case 3: Infinite Solutions (Dependent System)
**Input:**
```
3
1 2 3 6
2 4 6 12
3 6 9 18
```

**System of Equations:**
```
x₁ + 2x₂ + 3x₃ = 6
2x₁ + 4x₂ + 6x₃ = 12
3x₁ + 6x₂ + 9x₃ = 18
```

**Output:**
```
Infinite Solutions
```

**Explanation:** All three equations are scalar multiples of each other, representing the same plane in 3D space.

## 📂 Quick Links

| Resource | Description | Link |
|----------|-------------|------|
| 📄 **Source Code** | Complete C++ implementation | [gauss-elimination-method.cpp](gauss-elimination-method.cpp) |
| 📥 **Input File** | Sample test cases | [input.txt](input.txt) |
| 📤 **Output File** | Results with step-by-step solutions | [output.txt](output.txt) |

## 🚀 How to Run

1. **Compile the program:**
   ```bash
   g++ gauss-elimination-method.cpp -o gauss
   ```

2. **Prepare your input:**
   - Create or modify `input.txt` with your system(s) of equations
   - Format: First line = number of variables (n)
   - Next n lines = coefficients followed by the constant term

3. **Run the program:**
   ```bash
   ./gauss
   ```

4.  **View results:**
   - Check `output.txt` for detailed solutions and intermediate steps

## 📝 Input Format

```
n
a₁₁ a₁₂ ...  a₁ₙ b₁
a₂₁ a₂₂ ...  a₂ₙ b₂
...
aₙ₁ aₙ₂ ... aₙₙ bₙ
```

Where:
- `n` = number of variables/equations
- `aᵢⱼ` = coefficient of xⱼ in equation i
- `bᵢ` = constant term of equation i

## ✨ Features

- ✅ Partial pivoting for numerical stability
- ✅ Detects all three solution types (unique, no solution, infinite)
- ✅ Shows intermediate steps of elimination process
- ✅ Handles multiple test cases in a single run
- ✅ File-based input/output for easy testing
- ✅ Precision control (2 decimal places)

## 🔧 Configuration

Toggle intermediate step printing:
```cpp
bool printIntermediate = true; // Set to false to hide intermediate matrices
```

## 📚 References

- [Gaussian Elimination - Wikipedia](https://en.wikipedia.org/wiki/Gaussian_elimination)
- Numerical Methods For Engineers by Raymond Canale and Steven C. Chapra

---

**Part of the [Numerical Computing Suite](../) by [AbirHasanArko](https://github.com/AbirHasanArko)**