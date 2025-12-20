
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
    fout << fixed << setprecision(6);
    fout << "Number of points: " << n << endl;
    fout << "\nInput Points (x, y):" << endl;
    for(int i=0;i<n;i++) {
        fin >> x[i] >> y[i];
        fout << "(" << x[i] << ", " << y[i] << ")" << endl;
    }

    double sumX=0,sumX2=0,sumX3=0,sumX4=0,sumY=0,sumXY=0,sumX2Y=0;
    for(int i=0;i<n;i++){
        sumX += x[i];
        sumX2 += x[i]*x[i];
        sumX3 += x[i]*x[i]*x[i];
        sumX4 += x[i]*x[i]*x[i]*x[i];
        sumY += y[i];
        sumXY += x[i]*y[i];
        sumX2Y += x[i]*x[i]*y[i];
    }

    fout << "\nSum of x: " << sumX << endl;
    fout << "Sum of x^2: " << sumX2 << endl;
    fout << "Sum of x^3: " << sumX3 << endl;
    fout << "Sum of x^4: " << sumX4 << endl;
    fout << "Sum of y: " << sumY << endl;
    fout << "Sum of x*y: " << sumXY << endl;
    fout << "Sum of x^2*y: " << sumX2Y << endl;

    // Solve normal equations for quadratic fit: y = a + b*x + c*x^2
    // [ n    sumX   sumX2 ] [a]   [ sumY   ]
    // [sumX sumX2  sumX3 ] [b] = [ sumXY  ]
    // [sumX2 sumX3 sumX4] [c]   [ sumX2Y ]
    double A[3][4] = {
        {double(n), sumX, sumX2, sumY},
        {sumX, sumX2, sumX3, sumXY},
        {sumX2, sumX3, sumX4, sumX2Y}
    };

    // Gaussian elimination
    for(int i=0;i<3;i++){
        // Partial pivoting
        int maxRow = i;
        for(int k=i+1;k<3;k++){
            if(abs(A[k][i]) > abs(A[maxRow][i])) maxRow = k;
        }
        for(int k=i;k<4;k++) swap(A[maxRow][k], A[i][k]);
        // Eliminate
        for(int k=i+1;k<3;k++){
            double f = A[k][i]/A[i][i];
            for(int j=i;j<4;j++)
                A[k][j] -= f*A[i][j];
        }
    }
    // Back substitution
    double coeff[3];
    for(int i=2;i>=0;i--){
        coeff[i] = A[i][3];
        for(int j=i+1;j<3;j++)
            coeff[i] -= A[i][j]*coeff[j];
        coeff[i] /= A[i][i];
    }

    fout << "\nEquation of best fit quadratic polynomial (Least Squares):" << endl;
    fout << "y = " << coeff[0] << " + " << coeff[1] << "x + " << coeff[2] << "x^2" << endl;

    fin.close();
    fout.close();
    return 0;
}
