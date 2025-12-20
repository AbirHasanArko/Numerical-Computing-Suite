
#include <iostream>
#include <fstream>
#include <iomanip>
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
    double sumX=0, sumY=0, sumXY=0, sumX2=0;

    fout << fixed << setprecision(6);
    fout << "Number of points: " << n << endl;
    fout << "\nInput Points (x, y):" << endl;
    for(int i=0;i<n;i++){
        fin >> x[i] >> y[i];
        fout << "(" << x[i] << ", " << y[i] << ")" << endl;
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i]*y[i];
        sumX2 += x[i]*x[i];
    }

    fout << "\nSum of x: " << sumX << endl;
    fout << "Sum of y: " << sumY << endl;
    fout << "Sum of x*y: " << sumXY << endl;
    fout << "Sum of x^2: " << sumX2 << endl;

    double b = (n*sumXY - sumX*sumY) / (n*sumX2 - sumX*sumX);
    double a = (sumY - b*sumX) / n;

    fout << "\nEquation of best fit line (Least Squares):" << endl;
    fout << "y = " << a << " + " << b << "x" << endl;

    fin.close();
    fout.close();
    return 0;
}
