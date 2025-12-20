
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

int main() {
    int n;
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin || !fout) {
        cerr << "Error opening file." << endl;
        return 1;
    }

    fin >> n;
    double x[100], y[100];
    fout << fixed << setprecision(6);
    fout << "Number of points: " << n << endl;
    fout << "\nInput Points (x, y):" << endl;
    for(int i=0;i<n;i++) {
        fin >> x[i] >> y[i];
        fout << "(" << x[i] << ", " << y[i] << ")" << endl;
    }

    // Fit y = a * exp(bx) using linearization: ln(y) = ln(a) + b*x
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    double lnY[100];
    for(int i=0;i<n;i++) {
        lnY[i] = log(y[i]);
        sumX += x[i];
        sumY += lnY[i];
        sumXY += x[i]*lnY[i];
        sumX2 += x[i]*x[i];
    }

    fout << "\nSum of x: " << sumX << endl;
    fout << "Sum of ln(y): " << sumY << endl;
    fout << "Sum of x*ln(y): " << sumXY << endl;
    fout << "Sum of x^2: " << sumX2 << endl;

    double b = (n*sumXY - sumX*sumY) / (n*sumX2 - sumX*sumX);
    double A = (sumY - b*sumX) / n;
    double a = exp(A);

    fout << "\nEquation of best fit (Non-linear, y = a*exp(bx)):" << endl;
    fout << "y = " << a << " * exp(" << b << "x)" << endl;

    fin.close();
    fout.close();
    return 0;
}
