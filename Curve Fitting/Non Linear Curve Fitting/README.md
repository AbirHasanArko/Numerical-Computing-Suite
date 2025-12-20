## Non-Linear Curve Fitting

## 📑 Table of Contents
- [Introduction](#-introduction-2)
- [Common Models](#-common-models)
- [Linearization Technique](#-linearization-technique)
- [Algorithm](#-algorithm-2)
- [Applications](#-applications-2)

---

## 📖 Introduction
Non-linear curve fitting is used when data cannot be accurately
represented by linear or polynomial models.
Such equations are transformed into linear form
before applying the least squares method.

---

## 📌 Common Models
Some widely used non-linear models are:

- **Exponential curve**
y = a·e^(b·x)

- **Power curve**
y = a·x^b

---

## 🔄 Linearization Technique
For the exponential model:

y = a·e^(b·x)

Taking logarithm on both sides:

ln(y) = ln(a) + b·x

This converts the equation into a linear form suitable for least squares fitting.

---

## 🧾 Algorithm
1. Read the given data points
2. Apply logarithmic transformation
3. Convert the equation into linear form
4. Apply least squares method
5. Compute constants
6. Convert back to original non-linear equation
7. Display the fitted curve

---

## 🧪 Applications
- Population growth models
- Chemical reaction analysis
- Biological and economic data modeling