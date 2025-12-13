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

        // Gauss-Jordan Elimination
        int step = 1;
        for (int i = 0; i < n; i++)
        {
            // Partial Pivoting
            int maxRow = i;
            for (int k = i + 1; k < n; k++)
                if (fabs(a[k][i]) > fabs(a[maxRow][i]))
                    maxRow = k;
            
            if (i != maxRow)
            {
                swap(a[i], a[maxRow]);
                if (printIntermediate)
                {
                    fout << "\nAfter swapping row " << i + 1 << " with row " << maxRow + 1 << ":\n";
                    for (int r = 0; r < n; r++)
                    {
                        for (int c = 0; c <= n; c++)
                            fout << a[r][c] << "\t";
                        fout << "\n";
                    }
                }
            }

            if (fabs(a[i][i]) < 1e-12)
                continue;

            // Make diagonal element 1 (normalize pivot row)
            double pivot = a[i][i];
            for (int j = i; j <= n; j++)
                a[i][j] /= pivot;

            if (printIntermediate)
            {
                fout << "\nStep " << step++ << " - Making diagonal element a[" << i + 1 << "][" << i + 1 << "] = 1:\n";
                for (int r = 0; r < n; r++)
                {
                    for (int c = 0; c <= n; c++)
                        fout << a[r][c] << "\t";
                    fout << "\n";
                }
            }

            // Eliminate column i in ALL other rows (both above and below)
            for (int k = 0; k < n; k++)
            {
                if (k != i)
                {
                    double factor = a[k][i];
                    for (int j = i; j <= n; j++)
                        a[k][j] -= factor * a[i][j];
                }
            }

            // Print intermediate matrix if enabled
            if (printIntermediate)
            {
                fout << "\nStep " << step++ << " - Eliminating column " << i + 1 << " in all other rows:\n";
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
            if (! allZero)
                rank++;
        }

        if (! noSolution && rank < n)
            infiniteSolution = true;

        // Write output
        fout << "\n========================================\n";
        fout << "FINAL RESULT:\n";
        fout << "========================================\n";

        if (noSolution)
        {
            fout << "\nNo Solution\n";
            fout << "The system is inconsistent.\n";
        }
        else if (infiniteSolution)
        {
            fout << "\nInfinite Solutions\n";
            fout << "The system has dependent equations.\n";
        }
        else
        {
            fout << "\nUnique Solution\n";
            fout << "\nFinal Reduced Row Echelon Form (RREF):\n";
            for (int r = 0; r < n; r++)
            {
                for (int c = 0; c <= n; c++)
                    fout << a[r][c] << "\t";
                fout << "\n";
            }

            fout << "\nSolution:\n";
            for (int i = 0; i < n; i++)
            {
                x[i] = a[i][n];
                fout << "x" << i + 1 << " = " << x[i] << "\n";
            }
        }

        fout << "\n" << string(60, '=') << "\n\n";
    }

    fin.close();
    fout. close();
    cout << "All results written to output.txt\n";
    return 0;
}