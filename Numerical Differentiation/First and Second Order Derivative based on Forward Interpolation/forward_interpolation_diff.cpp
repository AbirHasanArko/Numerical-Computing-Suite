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

// Utility: Forward Difference Table
vector<vector<double>> forwardDiff(const vector<double> &y)
{
    int n = y.size();
    vector<vector<double>> d(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        d[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            d[i][j] = d[i + 1][j - 1] - d[i][j - 1];

    return d;
}

// Forward First Derivative
static double forwardFirstDerivative(const vector<double> &x,
                                     const vector<double> &y,
                                     double xp)
{
    double h = x[1] - x[0];
    auto d = forwardDiff(y);
    double u = (xp - x[0]) / h;
    double res = d[0][1];

    if (x.size() >= 3)
        res += ((2 * u - 1) / 2.0) * d[0][2];
    if (x.size() >= 4)
        res += ((3 * u * u - 6 * u + 2) / 6.0) * d[0][3];

    return res / h;
}

// Forward Second Derivative
static double forwardSecondDerivative(const vector<double> &x,
                                      const vector<double> &y,
                                      double xp)
{
    double h = x[1] - x[0];
    auto d = forwardDiff(y);
    double u = (xp - x[0]) / h;
    double res = d[0][2];

    if (x.size() >= 4)
        res += (u - 1) * d[0][3];
    if (x.size() >= 5)
        res += ((6 * u * u - 18 * u + 11) / 12.0) * d[0][4];
        
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
    fout << "Numerical Differentiation using Forward Interpolation\n";
    fout << "---------------------------------------------------\n";
    
    fout << "Given data points (x, y):\n";
    for (int i = 0; i < n; ++i)
        fout << "  (" << setw(8) << x[i] << ", " << setw(10) << y[i] << ")\n";
    fout << "\n";
    
    fout << "Point of differentiation (xp): " << xp << "\n\n";

    double first = forwardFirstDerivative(x, y, xp);
    double second = forwardSecondDerivative(x, y, xp);

    fout << "First Derivative at x = " << xp << " : " << first << "\n";
    fout << "Second Derivative at x = " << xp << " : " << second << "\n";
    fout << "\n";

    fin.close();
    fout.close();

    return 0;
}
