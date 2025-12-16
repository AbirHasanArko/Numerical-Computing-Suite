# Non-Linear Curve Fitting

---

## 📖 Introduction
Non-linear curve fitting is used when data cannot be represented by a
straight line or simple polynomial.
Common non-linear models include:

- Exponential curve: \( y = ae^{bx} \)
- Power curve: \( y = ax^b \)

These equations are transformed into linear form
using logarithmic transformations.

---

## 📌 Mathematical Transformation
Example (Exponential):

\[
y = ae^{bx}
\]

Taking log on both sides:

\[
\ln y = \ln a + bx
\]

This converts the equation into a linear form.

---

## 🧾 Algorithm Steps
1. Read the given data points
2. Transform the equation into linear form
3. Apply least squares method
4. Compute constants
5. Convert back to original equation
6. Display final curve equation

---

## ⚙️ Implementation Notes
- Requires mathematical transformation
- Accuracy depends on correct model selection

---

## 🧪 Usage Example
Refer to `input.txt` and `output.txt`.

---

## 📚 References
- S. S. Sastry  
- B. S. Grewal