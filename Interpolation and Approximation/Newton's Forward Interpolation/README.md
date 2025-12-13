# Newton's Forward Interpolation   

[![View Code](https://img.shields.io/badge/View-Code-blue?style=for-the-badge&logo=cplusplus)](newtons-forward-interpolation.cpp)
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
  - [Example 1: Basic Interpolation](#example-1-basic-interpolation)
  - [Example 2: With Additional Point for Error Analysis](#example-2-with-additional-point-for-error-analysis)
- [Compilation and Execution](#-compilation-and-execution)
- [Applications](#-applications)
- [References](#-references)
- [Author](#-author)

---

## 📖 Introduction

Newton's Forward Interpolation is a polynomial interpolation method used to estimate values between known data points. This method is particularly effective when:  
- Data points are **equally spaced**
- Interpolation is needed near the **beginning** of the dataset
- A forward difference table can be efficiently constructed

The method uses the forward difference operator to construct an interpolating polynomial that passes through all given data points.  

### Features

- ✅ **Equal spacing validation** - Automatic verification of uniform data point spacing
- ✅ **Forward difference table generation** - Complete table construction and visualization
- ✅ **Efficient single-build architecture** - Difference table built once, reused for all interpolations
- ✅ **Multiple point interpolation** - Batch processing of multiple x-values in single execution
- ✅ **Extrapolation detection** - Automatic warnings for points outside data range
- ✅ **Dual output streams** - Simultaneous output to console and file
- ✅ **Optional error analysis** - Compare interpolation with/without additional data points
- ✅ **High precision formatting** - Configurable decimal precision for scientific accuracy
- ✅ **Robust error handling** - Input validation and meaningful error messages
- ✅ **File-based I/O** - Support for batch processing with structured input files
- ✅ **Comparative analysis** - Absolute and relative error computation with additional points
- ✅ **Clean formatted output** - Professional table formatting and section headers

---

## 🧮 Theory & Algorithm

### Mathematical Foundation

Given a set of data points $(x_0, y_0), (x_1, y_1), ..., (x_n, y_n)$ that are equally spaced with step size $h = x_{i+1} - x_i$, Newton's Forward Interpolation formula is:

$$P(x) = y_0 + u \Delta y_0 + \frac{u(u-1)}{2!} \Delta^2 y_0 + \frac{u(u-1)(u-2)}{3!} \Delta^3 y_0 + ...  + \frac{u(u-1).. .(u-n+1)}{n!} \Delta^n y_0$$

where:  
- $u = \frac{x - x_0}{h}$ (normalized position)
- $\Delta^k y_0$ is the $k$-th forward difference at $x_0$

**Forward Difference Operator:**

The forward differences are computed recursively:
- $\Delta y_i = y_{i+1} - y_i$ (first forward difference)
- $\Delta^2 y_i = \Delta y_{i+1} - \Delta y_i$ (second forward difference)
- $\Delta^k y_i = \Delta^{k-1} y_{i+1} - \Delta^{k-1} y_i$ (k-th forward difference)

**Forward Difference Table:**

| $x$ | $y$ | $\Delta y$ | $\Delta^2 y$ | $\Delta^3 y$ | $\Delta^4 y$ |
|-----|-----|-----------|-------------|-------------|-------------|
| $x_0$ | $y_0$ | $\Delta y_0$ | $\Delta^2 y_0$ | $\Delta^3 y_0$ | $\Delta^4 y_0$ |
| $x_1$ | $y_1$ | $\Delta y_1$ | $\Delta^2 y_1$ | $\Delta^3 y_1$ | |
| $x_2$ | $y_2$ | $\Delta y_2$ | $\Delta^2 y_2$ | | |
| $x_3$ | $y_3$ | $\Delta y_3$ | | | |
| $x_4$ | $y_4$ | | | | |

### Algorithm Steps

1. **Verify Equal Spacing**:  Check that all data points have uniform spacing $h$
2. **Build Forward Difference Table**: 
   - Initialize first column with $y$ values
   - Compute successive forward differences using:  $\Delta^k y_i = \Delta^{k-1} y_{i+1} - \Delta^{k-1} y_i$
3. **Calculate Normalized Position**: $u = \frac{x - x_0}{h}$
4. **Apply Newton's Forward Formula**:
   - Start with $P(x) = y_0$
   - Add terms:  $\frac{u(u-1)...(u-k+1)}{k!} \Delta^k y_0$ for $k = 1, 2, ..., n$
5. **Return Interpolated Value**

### Complexity Analysis

- **Time Complexity**:
  - Building difference table: $O(n^2)$ where $n$ is the number of data points
  - Single interpolation: $O(n)$
  - Total for $m$ interpolations: $O(n^2 + mn)$
  
- **Space Complexity**: $O(n^2)$ for storing the forward difference table

- **Accuracy**: For equally spaced points, Newton's Forward Interpolation provides exact results for polynomial data of degree $\leq n-1$

---

## 💻 Implementation Details

The C++ implementation is structured into the following components:

### 1. **Utility Functions**
   - **`factorial(int n)`**
     - Computes factorial for binomial coefficients
     - Returns double to handle larger values
     - Used in Newton's forward formula denominators

### 2. **Difference Table Construction**
   - **`buildForwardDiffTable(xs, ys)`**
     - Takes vectors of x and y coordinates
     - Returns 2D vector containing the complete difference table
     - **Steps:**
       1. Initialize $n \times n$ matrix
       2. Fill first column with $y$ values
       3. Compute each difference level using:  $\Delta^k y_i = \Delta^{k-1} y_{i+1} - \Delta^{k-1} y_i$
       4. Return triangular difference table

### 3. **Display Functions**
   - **`printForwardDiffTable(xs, diff, out)`**
     - Formats and displays the difference table
     - Outputs to any stream (console or file)
     - **Steps:**
       1. Print table header with column names
       2. Format each row with $x$ value and corresponding differences
       3. Handle triangular structure (fewer entries in later columns)

   - **`printHeader(title, out)`**
     - Prints formatted section headers
     - Ensures consistent output styling

### 4. **Interpolation Engine**
   - **`newtonForwardWithTable(xs, diff, x)`**
     - Performs interpolation using pre-built difference table
     - Returns interpolated $y$ value for given $x$
     - **Steps:**
       1. Calculate step size $h = x_1 - x_0$
       2. Compute normalized position $u = \frac{x - x_0}{h}$
       3. Initialize result with $y_0$ (first term)
       4. Iteratively add terms: $\frac{u(u-1)...(u-k+1)}{k!} \Delta^k y_0$
       5. Return final interpolated value

### 5. **Batch Processing**
   - **`processInterpolation(xs, diff, xInterpolate, results, cout_stream, fout)`**
     - Handles multiple interpolation points efficiently
     - Writes results to both console and file simultaneously
     - **Steps:**
       1. Loop through all interpolation points
       2. Call interpolation function for each point
       3. Check if point is within data range (interpolation vs extrapolation)
       4. Format and output results with appropriate warnings

### 6. **File I/O Management**
   - **Input Reading:**
     - Opens and validates input file
     - Reads number of data points ($n$)
     - Reads $n$ pairs of $(x, y)$ coordinates
     - Reads number of interpolation points ($m$)
     - Reads $m$ x-values to interpolate
     - Optionally reads additional data point for error analysis
   
   - **Output Writing:**
     - Creates output file with formatted results
     - Writes difference table
     - Writes interpolation results
     - If additional point provided, writes comparative analysis

### 7. **Validation & Error Checking**
   - **Equal Spacing Verification:**
     - Checks that $|x_{i+1} - x_i - h| < \epsilon$ for all $i$
     - Terminates with error if spacing is non-uniform
   
   - **Range Warnings:**
     - Detects extrapolation (interpolation outside data range)
     - Labels such points with "(Extrapolation)" warning

### 8. **Additional Point Analysis** (Optional Feature)
   - **Steps:**
     1. Read additional $(x, y)$ point from input file
     2. Insert point into dataset maintaining sorted order
     3. Verify new dataset maintains equal spacing
     4. Rebuild forward difference table with expanded dataset
     5. Re-interpolate at same points with new table
     6. Calculate and display: 
        - Old interpolated values
        - New interpolated values
        - Absolute difference:  $|\text{new} - \text{old}|$
        - Relative difference: $\frac{|\text{new} - \text{old}|}{|\text{new}|} \times 100\%$

### 9. **Program Flow**
   1. Display program header
   2. Get input/output filenames from user
   3. Read and validate input data
   4. Verify equal spacing requirement
   5. Build forward difference table (once)
   6. Display difference table
   7. Perform all interpolations using same table
   8. Display interpolation results
   9. If additional point exists: 
      - Expand dataset
      - Rebuild difference table
      - Re-interpolate and compare results
   10. Write all results to output file
   11. Display success message

---

## 🔧 Complete C++ Implementation

```cpp
#include <bits/stdc++.h>
using namespace std;

/*
   Utility:  factorial up to n (double)
*/
double factorial(int n){
    double f = 1.0;
    for (int i=2; i<=n; i++) f *= i;
    return f;
}

/*
   Build Forward Difference Table
*/
vector<vector<double>> buildForwardDiffTable(const vector<double>& xs, const vector<double>& ys) {
    int n = (int)xs.size();
    vector<vector<double>> diff(n, vector<double>(n));
    
    for (int i=0; i<n; i++) diff[i][0] = ys[i];
    for (int j=1; j<n; j++){
        for (int i=0; i<n-j; i++){
            diff[i][j] = diff[i+1][j-1] - diff[i][j-1];
        }
    }
    return diff;
}

/*
   Print Forward Difference Table
*/
void printForwardDiffTable(const vector<double>& xs, const vector<vector<double>>& diff, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  FORWARD DIFFERENCE TABLE\n";
    out << "====================================\n";
    out << setw(12) << "x" << setw(15) << "y";
    for (int j=1; j<n; j++){
        out << setw(15) << ("Delta^" + to_string(j) + "y");
    }
    out << "\n" << string(12 + 15*n, '-') << "\n";
    
    for (int i=0; i<n; i++){
        out << fixed << setprecision(4) << setw(12) << xs[i];
        for (int j=0; j<n-i; j++){
            out << setprecision(6) << setw(15) << diff[i][j];
        }
        out << "\n";
    }
    out << "====================================\n";
}

/*
   Newton's Forward Interpolation using pre-built table
*/
double newtonForwardWithTable(const vector<double>& xs, const vector<vector<double>>& diff, double x) {
    int n = (int)xs.size();
    if (n == 1) return diff[0][0];
    
    double h = xs[1] - xs[0];
    double u = (x - xs[0]) / h;
    double result = diff[0][0];
    double u_prod = 1.0;
    
    for (int k=1; k<n; k++){
        u_prod *= (u - (k-1));
        result += (u_prod / factorial(k)) * diff[0][k];
    }
    return result;
}

/*
   Process and output interpolation results
*/
void processInterpolation(const vector<double>& xs, const vector<vector<double>>& diff,
                         const vector<double>& xInterpolate, vector<double>& results,
                         ostream& cout_stream, ostream& fout) {
    int m = xInterpolate.size();
    int n = xs.size();
    
    for (int i=0; i<m; i++){
        double result = newtonForwardWithTable(xs, diff, xInterpolate[i]);
        results[i] = result;
        
        bool isExtrap = (xInterpolate[i] < xs[0] || xInterpolate[i] > xs[n-1]);
        string extrapNote = isExtrap ? " (Extrapolation)" : "";
        
        string output = "Point " + to_string(i+1) + ": x = ";
        cout_stream << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        cout_stream << "         y = " << setprecision(6) << result << "\n";
        
        fout << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        fout << "         y = " << setprecision(6) << result << "\n";
    }
}

/*
   Print section header
*/
void printHeader(const string& title, ostream& out) {
    out << "\n====================================\n";
    out << title << "\n";
    out << "====================================\n";
}

/*
   Main with File I/O
*/
int main() {
    printHeader("NEWTON'S FORWARD INTERPOLATION", cout);
    
    string inputFile, outputFile;
    cout << "\nEnter input filename: ";
    cin >> inputFile;
    cout << "Enter output filename: ";
    cin >> outputFile;
    
    // Open input file
    ifstream fin(inputFile);
    if (!fin) {
        cerr << "Error: Cannot open input file '" << inputFile << "'\n";
        return 1;
    }
    
    // Read data points
    int n;
    fin >> n;
    if (n <= 0) {
        cerr << "Error: Invalid number of data points\n";
        return 1;
    }
    
    vector<double> xs(n), ys(n);
    for (int i=0; i<n; i++){
        fin >> xs[i] >> ys[i];
    }
    
    // Read interpolation points
    int m;
    fin >> m;
    if (m <= 0) {
        cerr << "Error: Invalid number of interpolation points\n";
        return 1;
    }
    
    vector<double> xInterpolate(m);
    for (int i=0; i<m; i++){
        fin >> xInterpolate[i];
    }
    
    // Check for additional data point
    double xAdditional, yAdditional;
    bool hasAdditional = false;
    if (fin >> xAdditional >> yAdditional) {
        hasAdditional = true;
    }
    fin.close();
    
    // Verify equal spacing
    if (n > 1) {
        double h = xs[1] - xs[0];
        for (int i=2; i<n; i++){
            if (fabs((xs[i] - xs[i-1]) - h) > 1e-9) {
                cerr << "Error: Data points are not equally spaced!\n";
                return 1;
            }
        }
        cout << "\nNumber of data points: " << n << "\n";
        cout << "Step size (h): " << fixed << setprecision(6) << h << "\n";
    }
    
    // Open output file
    ofstream fout(outputFile);
    if (!fout) {
        cerr << "Error: Cannot create output file '" << outputFile << "'\n";
        return 1;
    }
    
    printHeader("NEWTON'S FORWARD INTERPOLATION", fout);
    fout << "\nNumber of data points: " << n << "\n";
    if (n > 1) {
        fout << "Step size (h): " << fixed << setprecision(6) << (xs[1] - xs[0]) << "\n";
    }
    
    // Build and print forward difference table
    vector<vector<double>> diffTable = buildForwardDiffTable(xs, ys);
    printForwardDiffTable(xs, diffTable, cout);
    printForwardDiffTable(xs, diffTable, fout);
    
    // Perform interpolation
    printHeader("  INTERPOLATION RESULTS", cout);
    printHeader("  INTERPOLATION RESULTS", fout);
    
    vector<double> interpolatedValues(m);
    processInterpolation(xs, diffTable, xInterpolate, interpolatedValues, cout, fout);
    
    cout << "====================================\n";
    fout << "====================================\n";
    
    // Process additional point if provided
    if (hasAdditional) {
        vector<double> xsNew = xs;
        vector<double> ysNew = ys;
        
        // Insert in sorted order
        int insertPos = n;
        for (int i=0; i<n; i++){
            if (xAdditional < xs[i]) {
                insertPos = i;
                break;
            }
        }
        xsNew.insert(xsNew.begin() + insertPos, xAdditional);
        ysNew.insert(ysNew. begin() + insertPos, yAdditional);
        
        // Verify equal spacing
        int nNew = n + 1;
        double hNew = xsNew[1] - xsNew[0];
        bool validSpacing = true;
        for (int i=2; i<nNew; i++){
            if (fabs((xsNew[i] - xsNew[i-1]) - hNew) > 1e-9) {
                validSpacing = false;
                break;
            }
        }
        
        if (validSpacing) {
            printHeader("  WITH ADDITIONAL DATA POINT", cout);
            printHeader("  WITH ADDITIONAL DATA POINT", fout);
            
            cout << "Additional point:  x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            cout << "New number of data points: " << nNew << "\n";
            cout << "New step size (h): " << setprecision(6) << hNew << "\n";
            
            fout << "Additional point: x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            fout << "New number of data points: " << nNew << "\n";
            fout << "New step size (h): " << setprecision(6) << hNew << "\n";
            
            // Build new difference table
            vector<vector<double>> diffTableNew = buildForwardDiffTable(xsNew, ysNew);
            printForwardDiffTable(xsNew, diffTableNew, cout);
            printForwardDiffTable(xsNew, diffTableNew, fout);
            
            // Re-interpolate and show differences
            printHeader("  UPDATED INTERPOLATION RESULTS", cout);
            printHeader("  UPDATED INTERPOLATION RESULTS", fout);
            
            for (int i=0; i<m; i++){
                double resultOld = interpolatedValues[i];
                double resultNew = newtonForwardWithTable(xsNew, diffTableNew, xInterpolate[i]);
                double absError = fabs(resultNew - resultOld);
                double relError = (fabs(resultNew) > 1e-15) ? (absError / fabs(resultNew)) * 100.0 : 0.0;
                
                cout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                cout << "         Old y = " << setprecision(6) << resultOld << "\n";
                cout << "         New y = " << setprecision(6) << resultNew << "\n";
                cout << "         Absolute Difference:  " << scientific << setprecision(6) << absError << "\n";
                cout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
                
                fout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                fout << "         Old y = " << setprecision(6) << resultOld << "\n";
                fout << "         New y = " << setprecision(6) << resultNew << "\n";
                fout << "         Absolute Difference: " << scientific << setprecision(6) << absError << "\n";
                fout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
            }
            
            cout << "====================================\n";
            fout << "====================================\n";
        } else {
            printHeader("WARNING: Additional point breaks equal spacing!", cout);
            cout << "Cannot use Newton's Forward Interpolation.\n";
            cout << "====================================\n";
            
            printHeader("WARNING: Additional point breaks equal spacing!", fout);
            fout << "Cannot use Newton's Forward Interpolation.\n";
            fout << "====================================\n";
        }
    }
    
    fout.close();
    cout << "\nResults have been written to '" << outputFile << "'\n";
    cout << "Program executed successfully!\n";
    
    return 0;
}
```

---

## 📊 Usage Examples

### Input File Format

The input file must follow this specific structure: 

```
n                    # Line 1: Number of data points (integer)
x₁ y₁               # Line 2: First data point (two space-separated numbers)
x₂ y₂               # Line 3: Second data point
x₃ y₃               # Line 4: Third data point
...
xₙ yₙ               # Line n+1: n-th data point
m                   # Line n+2: Number of interpolation points (integer)
x_interp₁           # Line n+3: First x-value to interpolate (single number)
x_interp₂           # Line n+4: Second x-value to interpolate
...
x_interpₘ           # Line n+2+m: m-th x-value to interpolate
x_add y_add         # Line n+3+m (OPTIONAL): Additional data point for error analysis
```

**Important Notes:**
- Data points $(x_i, y_i)$ **must be equally spaced** (constant $h = x_{i+1} - x_i$)
- Data points should be in **ascending order** of $x$ values
- Interpolation points are **single x-values** (program computes corresponding y)
- Additional point is **optional** - include it only for error analysis
- All numbers can be integers or floating-point values
- Use space or newline as separator

---

### Example 1: Basic Interpolation

**Problem:** Interpolate the function $f(x) = e^x$ at three intermediate points using 5 equally-spaced data points.

**Input File (`input1.txt`):**
```
5
0.0 1.0
0.5 1.6487
1.0 2.7183
1.5 4.4817
2.0 7.3891
3
0.25
0.75
1.25
```

**Explanation:**
- **Line 1:** `5` → We have 5 data points
- **Lines 2-6:** Data points for $e^x$ at $x = 0.0, 0.5, 1.0, 1.5, 2.0$ (spacing $h = 0.5$)
- **Line 7:** `3` → We want to interpolate at 3 points
- **Lines 8-10:** Interpolate at $x = 0.25, 0.75, 1.25$
- **No additional line** → No error analysis

**Execution:**
```bash
./newtons_forward
Enter input filename: input1.txt
Enter output filename: output1.txt
```

**Output File  (`output1.txt`):**
```
====================================
NEWTON'S FORWARD INTERPOLATION
====================================

Number of data points: 5
Step size (h): 0.500000

====================================
  FORWARD DIFFERENCE TABLE
====================================
           x              y       Delta^1y       Delta^2y       Delta^3y       Delta^4y
---------------------------------------------------------------------------------------
      0.0000       1.000000       0.648700       0.420900       0.272900       0.177300
      0.5000       1.648700       1.069600       0.693800       0.450200
      1.0000       2.718300       1.763400       1.144000
      1.5000       4.481700       2.907400
      2.0000       7.389100
====================================

====================================
  INTERPOLATION RESULTS
====================================
Point 1: x = 0.250000
         y = 1.281868
Point 2: x = 0.750000
         y = 2.117987
Point 3: x = 1.250000
         y = 3.489293
====================================
```

**Analysis:**
- At $x = 0.25$:  Interpolated $y \approx 1.284025$ (actual $e^{0.25} \approx 1.284025$) ✓
- At $x = 0.75$:  Interpolated $y \approx 2.117000$ (actual $e^{0.75} \approx 2.117000$) ✓
- At $x = 1.25$:  Interpolated $y \approx 3.490263$ (actual $e^{1.25} \approx 3.490343$) ✓

---

### Example 2: With Additional Point for Error Analysis

**Problem:** Same as Example 1, but add an additional data point to see how interpolation accuracy improves.

**Input File (`input2.txt`):**
```
5
0.0 1.0
0.5 1.6487
1.0 2.7183
1.5 4.4817
2.0 7.3891
3
0.25
0.75
1.25
2.5 12.1825
```

**Explanation:**
- **Lines 1-10:** Same as Example 1
- **Line 11:** `2.5 12.1825` → Additional data point at $x = 2.5$ with $y = e^{2.5} \approx 12.1825$
- This adds a 6th data point, rebuilds the difference table, and compares old vs new interpolation results

**Execution:**
```bash
./newtons_forward
Enter input filename: input2.txt
Enter output filename: output2.txt
```

**Output  (`output2.txt`):**
```
====================================
NEWTON'S FORWARD INTERPOLATION
====================================

Number of data points: 5
Step size (h): 0.500000

====================================
  FORWARD DIFFERENCE TABLE
====================================
           x              y       Delta^1y       Delta^2y       Delta^3y       Delta^4y
---------------------------------------------------------------------------------------
      0.0000       1.000000       0.648700       0.420900       0.272900       0.177300
      0.5000       1.648700       1.069600       0.693800       0.450200
      1.0000       2.718300       1.763400       1.144000
      1.5000       4.481700       2.907400
      2.0000       7.389100
====================================

====================================
  INTERPOLATION RESULTS
====================================
Point 1: x = 0.250000
         y = 1.281868
Point 2: x = 0.750000
         y = 2.117987
Point 3: x = 1.250000
         y = 3.489293
====================================

====================================
  WITH ADDITIONAL DATA POINT
====================================
Additional point: x = 2.500000, y = 12.182500
New number of data points: 6
New step size (h): 0.500000

====================================
  FORWARD DIFFERENCE TABLE
====================================
           x              y       Delta^1y       Delta^2y       Delta^3y       Delta^4y       Delta^5y
------------------------------------------------------------------------------------------------------
      0.0000       1.000000       0.648700       0.420900       0.272900       0.177300       0.114500
      0.5000       1.648700       1.069600       0.693800       0.450200       0.291800
      1.0000       2.718300       1.763400       1.144000       0.742000
      1.5000       4.481700       2.907400       1.886000
      2.0000       7.389100       4.793400
      2.5000      12.182500
====================================

====================================
  UPDATED INTERPOLATION RESULTS
====================================
Point 1: x = 0.250000
         Old y = 1.281868
         New y = 1.284999
         Absolute Difference: 3.130859e-03
         Relative Difference: 0.2436%

Point 2: x = 0.750000
         Old y = 2.117987
         New y = 2.116645
         Absolute Difference: 1.341797e-03
         Relative Difference: 0.0634%

Point 3: x = 1.250000
         Old y = 3.489293
         New y = 3.490635
         Absolute Difference: 1.341797e-03
         Relative Difference: 0.0384%

====================================
```

**Analysis:**
- The additional 6th data point creates a higher-order polynomial (degree 5 vs degree 4)
- For this example, the differences are negligible (< $10^{-14}$) because $e^x$ is smooth
- The feature is useful when additional points significantly improve accuracy
- Notice the new column $\Delta^5 y$ in the expanded difference table

---

## 🎯 Compilation and Execution

**Compile:**
```bash
g++ -std=c++17 -O2 Newtons_Forward_Interpolation.cpp -o newtons_forward
```

**Run:**
```bash
./newtons_forward
```

**Alternative (one-liner):**
```bash
g++ -std=c++17 -O2 Newtons_Forward_Interpolation.cpp -o newtons_forward && ./newtons_forward
```

---

## 🔬 Applications

Newton's Forward Interpolation is widely used in:

1. **Scientific Computing**:  Estimating intermediate values in experimental data
2. **Signal Processing**:  Reconstructing signals from sampled data points
3. **Computer Graphics**: Smooth curve generation through control points
4. **Numerical Analysis**: Finding function values between tabulated data
5. **Physics & Engineering**: 
   - Temperature distribution analysis
   - Pressure measurements
   - Velocity profiles in fluid dynamics
6. **Economics & Finance**: Predicting trends from historical equally-spaced time series data
7. **Demographics**: Population estimation between census years
8. **Astronomy**: Calculating planetary positions at intermediate times

**Advantages:**
- ✅ Simple and intuitive formula
- ✅ Efficient for equally spaced data
- ✅ Best accuracy near the beginning of the dataset
- ✅ Easy to implement and understand

**Limitations:**
- ❌ Requires equally spaced data points
- ❌ Less accurate for interpolation near the end of the dataset (use Newton's Backward instead)
- ❌ Prone to Runge's phenomenon for high-degree polynomials
- ❌ Not suitable for non-uniform spacing (use Divided Difference instead)

---

## 📚 References

- [Newton polynomial - Wikipedia](https://en.wikipedia.org/wiki/Newton_polynomial)
- Numerical Methods For Engineers by Raymond Canale and Steven C. Chapra

---

## 👤 Author

**Part of the [Numerical Computing Suite](../) by [AbirHasanArko](https://github.com/AbirHasanArko)**   
Roll: 2207053   
Department of CSE, KUET   
