Simpson's Three-Eighths (3/8) Rule

[View Code] (simpsons-three-eighths-rule.cpp)
[View Input] (input.txt)
[View Output] (output.txt)

--------------------------------------------------

Table of Contents
1. Introduction
2. Mathematical Formula
3. Validity Condition
4. Algorithm Steps
5. Implementation Notes
6. Usage Example
7. References

--------------------------------------------------

Introduction
Simpson’s Three-Eighths (3/8) Rule approximates the definite integral by fitting cubic (3rd-degree)
polynomials over groups of three consecutive subintervals. It is especially useful when the
number of subintervals is a multiple of 3 and provides good accuracy for smooth functions.

--------------------------------------------------

Mathematical Formula
Divide the interval [a, b] into n equal subintervals (n must be a multiple of 3).

Step size:
h = (b - a) / n

Points:
xi = a + i·h , where i = 0, 1, 2, …, n

Formula:
∫ab f(x) dx ≈ (3h / 8) [ y0 + yn
+ 3(y1 + y2 + y4 + y5 + …)
+ 2(y3 + y6 + y9 + …) ]

where yi = f(xi)

--------------------------------------------------

Validity Condition
- n must be a multiple of 3
- Data points must be equally spaced

--------------------------------------------------

Algorithm Steps
1. Read n
2. Read a and b
3. Read y0 to yn (total n+1 values)
4. Check if n is divisible by 3
5. Compute h = (b - a) / n
6. Compute:
   - sum_3 = sum of yi where i is NOT divisible by 3
   - sum_2 = sum of yi where i IS divisible by 3 (excluding y0 and yn)
7. Apply Simpson’s 3/8 formula
8. Display the approximate integral value

--------------------------------------------------

Implementation Notes
- Program supports multiple test cases
- Output is formatted to match numerical methods lab requirements

--------------------------------------------------

Usage Example
See input.txt and output.txt

--------------------------------------------------

References
1. S. S. Sastry – Introductory Methods of Numerical Analysis
2. B. S. Grewal – Numerical Methods in Engineering and Science
