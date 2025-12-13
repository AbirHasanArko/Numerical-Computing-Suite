#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin)
    {
        cerr << "Error: input.txt not found!\n";
        return 1;
    }

    fout << fixed << setprecision(2);

    bool printIntermediate = true; // toggle intermediate steps

    int n;
    while (fin >> n)
    {
        vector<vector<double>> a(n, vector<double>(n + 1));
        vector<double> x(n);

        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                fin >> a[i][j];

        // Print the original system
        fout << "\nInput system:\n";
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (j > 0 && a[i][j] >= 0) fout << "+";
                fout << a[i][j] << "x" << j + 1 << " ";
            }
            fout << "= " << a[i][n] << "\n";
        }

        // Forward Elimination
        for (int i = 0; i < n - 1; i++)
        {
            int maxRow = i;
            for (int k = i + 1; k < n; k++)
                if (fabs(a[k][i]) > fabs(a[maxRow][i]))
                    maxRow = k;
            swap(a[i], a[maxRow]);

            if (fabs(a[i][i]) < 1e-12)
                continue;

            for (int k = i + 1; k < n; k++)
            {
                double factor = a[k][i] / a[i][i];
                for (int j = i; j <= n; j++)
                    a[k][j] -= factor * a[i][j];
            }

            // Print intermediate matrix if enabled
            if (printIntermediate)
            {
                fout << "\nAfter step " << i + 1 << ":\n";
                for (int r = 0; r < n; r++)
                {
                    for (int c = 0; c <= n; c++)
                        fout << a[r][c] << "\t";
                    fout << "\n";
                }
            }
        }

        // Detect solution type
        bool noSolution = false, infiniteSolution = false;
        int rank = 0;

        for (int i = 0; i < n; i++)
        {
            bool allZero = true;
            for (int j = 0; j < n; j++)
                if (fabs(a[i][j]) > 1e-12) allZero = false;

            if (allZero && fabs(a[i][n]) > 1e-12)
            {
                noSolution = true;
                break;
            }
            if (!allZero)
                rank++;
        }

        if (!noSolution && rank < n)
            infiniteSolution = true;

        // Write output
        if (noSolution)
            fout << "\nNo Solution\n";
        else if (infiniteSolution)
            fout << "\nInfinite Solutions\n";
        else
        {
            fout << "\nUnique Solution\n";
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = a[i][n];
                for (int j = i + 1; j < n; j++)
                    x[i] -= a[i][j] * x[j];
                x[i] /= a[i][i];
            }

            fout << "Solution:\n";
            for (int i = 0; i < n; i++)
                fout << "x" << i + 1 << " = " << x[i] << "\n";
        }

        fout << "\n-------------------------------------------\n\n";
    }

    fin.close();
    fout.close();
    cout << "All results written to output.txt\n";
    return 0;
}