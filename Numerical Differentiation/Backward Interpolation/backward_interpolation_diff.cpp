#include <bits/stdc++.h>
using namespace std;

// Utility: factorial
static double factorial(int n)
{
    double f = 1.0;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

// Utility: Backward Difference Table
vector<vector<double>> backwardDiff(const vector<double> &y)
{
    int n = y.size();

    vector<vector<double>> b(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        b[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = n - 1; i >= j; i--)
            b[i][j] = b[i][j - 1] - b[i - 1][j - 1];

    return b;
}

// Backward First Derivative
static double backwardFirstDerivative(const vector<double> &x,
                                      const vector<double> &y,
                                      double xp)
{
    int n = x.size();

    double h = x[1] - x[0];
    auto b = backwardDiff(y);
    double v = (xp - x[n - 1]) / h;
    double res = b[n - 1][1];

    if (n >= 3)
        res += ((2 * v + 1) / 2.0) * b[n - 1][2];
    if (n >= 4)
        res += ((3 * v * v + 6 * v + 2) / 6.0) * b[n - 1][3];

    return res / h;
}

// Backward Second Derivative
static double backwardSecondDerivative(const vector<double> &x,
                                       const vector<double> &y,
                                       double xp)
{
    int n = x.size();

    double h = x[1] - x[0];
    auto b = backwardDiff(y);
    double v = (xp - x[n - 1]) / h;
    double res = b[n - 1][2];

    if (n >= 4)
        res += (v + 1) * b[n - 1][3];
    if (n >= 5)
        res += ((6 * v * v + 18 * v + 11) / 12.0) * b[n - 1][4];
        
    return res / (h * h);
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if (!fin.is_open())
    {
        cerr << "Error: Could not open input.txt" << endl;
        return 1;
    }
    if (!fout.is_open())
    {
        cerr << "Error: Could not open output.txt" << endl;
        return 1;
    }

    int n;
    fin >> n;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; ++i)
        fin >> x[i];
    for (int i = 0; i < n; ++i)
        fin >> y[i];

    double xp;
    fin >> xp;

    fout << fixed << setprecision(6);
    fout << "Numerical Differentiation using Backward Interpolation\n";
    fout << "---------------------------------------------------\n";
    
    fout << "Given data points (x, y):\n";
    for (int i = 0; i < n; ++i)
        fout << "  (" << setw(8) << x[i] << ", " << setw(10) << y[i] << ")\n";
    fout << "\n";

    fout << "Point of differentiation (xp): " << xp << "\n\n";

    double first = backwardFirstDerivative(x, y, xp);
    double second = backwardSecondDerivative(x, y, xp);
    fout << "First Derivative at x = " << xp << " : " << first << "\n";
    fout << "Second Derivative at x = " << xp << " : " << second << "\n";
    fout << "\n";

    fin.close();
    fout.close();
    
    return 0;
}
