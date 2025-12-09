#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin)
    {
        cerr << "Error:  input.txt not found!\n";
        return 1;
    }

    fout << fixed << setprecision(4);

    bool printIntermediate = true; // toggle intermediate steps

    int n;
    while (fin >> n)
    {
        vector<vector<double>> aug(n, vector<double>(n + 1));
        
        // Read augmented matrix
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                fin >> aug[i][j];

        // Print the original system
        fout << "\n========================================\n";
        fout << "Input system:\n";
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (j > 0 && aug[i][j] >= 0) fout << "+";
                fout << aug[i][j] << "x" << j + 1 << " ";
            }
            fout << "= " << aug[i][n] << "\n";
        }
        fout << "========================================\n";

        // Separate A and b from augmented matrix
        vector<vector<double>> A(n, vector<double>(n));
        vector<double> b(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                A[i][j] = aug[i][j];
            b[i] = aug[i][n];
        }

        // Initialize L and U matrices
        vector<vector<double>> L(n, vector<double>(n, 0));
        vector<vector<double>> U(n, vector<double>(n, 0));

        fout << "\nPerforming LU Decomposition...\n";

        bool decompositionFailed = false;
        int failedAtStep = -1;

        // LU Decomposition using Doolittle's method
        for (int i = 0; i < n; i++)
        {
            // Upper Triangular Matrix U
            for (int k = i; k < n; k++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += L[i][j] * U[j][k];
                U[i][k] = A[i][k] - sum;
            }

            // Check for zero pivot
            if (fabs(U[i][i]) < 1e-12)
            {
                decompositionFailed = true;
                failedAtStep = i;
                break;
            }

            // Lower Triangular Matrix L
            for (int k = i; k < n; k++)
            {
                if (i == k)
                    L[i][i] = 1; // diagonal element is 1
                else
                {
                    double sum = 0;
                    for (int j = 0; j < i; j++)
                        sum += L[k][j] * U[j][i];
                    L[k][i] = (A[k][i] - sum) / U[i][i];
                }
            }

            // Print intermediate matrices if enabled
            if (printIntermediate)
            {
                fout << "\nAfter step " << i + 1 << ":\n";
                fout << "L matrix:\n";
                for (int r = 0; r < n; r++)
                {
                    for (int c = 0; c < n; c++)
                        fout << setw(10) << L[r][c] << " ";
                    fout << "\n";
                }
                fout << "U matrix:\n";
                for (int r = 0; r < n; r++)
                {
                    for (int c = 0; c < n; c++)
                        fout << setw(10) << U[r][c] << " ";
                    fout << "\n";
                }
                fout << "---------------------------------------------\n";
            }
        }

        // Calculate determinant (product of diagonal of U)
        double detU = 1;
        for (int i = 0; i < n; i++)
            detU *= U[i][i];

        // Detect solution type
        bool noSolution = false;
        bool infiniteSolution = false;

        if (decompositionFailed || fabs(detU) < 1e-12)
        {
            // Matrix is singular - check for no solution vs infinite solutions
            // Continue decomposition as much as possible to check consistency
            
            // Perform forward substitution to check consistency
            vector<double> y(n, 0);
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += L[i][j] * y[j];
                
                double expectedY = b[i] - sum;
                
                // Check if this row is all zeros in U
                bool rowIsZero = true;
                for (int j = i; j < n; j++)
                {
                    if (fabs(U[i][j]) > 1e-12)
                    {
                        rowIsZero = false;
                        break;
                    }
                }
                
                if (rowIsZero && fabs(expectedY) > 1e-12)
                {
                    noSolution = true;
                    break;
                }
                
                if (!rowIsZero && fabs(U[i][i]) > 1e-12)
                    y[i] = expectedY;
            }

            if (! noSolution)
                infiniteSolution = true;
        }

        // Output results
        fout << "\n========================================\n";
        fout << "FINAL RESULT:\n";
        fout << "========================================\n";

        if (noSolution)
        {
            fout << "\nNo Solution\n";
            fout << "The system is inconsistent.\n";
            fout << "Determinant of U = " << detU << " (approximately 0)\n";
        }
        else if (infiniteSolution)
        {
            fout << "\nInfinite Solutions\n";
            fout << "The system has dependent equations.\n";
            fout << "Determinant of U = " << detU << " (approximately 0)\n";
        }
        else
        {
            fout << "\nUnique Solution\n";
            fout << "Determinant of U = " << detU << "\n\n";

            // Print final L and U matrices
            fout << "Final L matrix (Lower Triangular):\n";
            for (int r = 0; r < n; r++)
            {
                for (int c = 0; c < n; c++)
                    fout << setw(10) << L[r][c] << " ";
                fout << "\n";
            }

            fout << "\nFinal U matrix (Upper Triangular):\n";
            for (int r = 0; r < n; r++)
            {
                for (int c = 0; c < n; c++)
                    fout << setw(10) << U[r][c] << " ";
                fout << "\n";
            }

            // Forward substitution to solve L*y = b
            vector<double> y(n);
            fout << "\n--- Forward Substitution (L*y = b) ---\n";
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += L[i][j] * y[j];
                y[i] = b[i] - sum;
                fout << "y" << i + 1 << " = " << y[i] << "\n";
            }

            // Back substitution to solve U*x = y
            vector<double> x(n);
            fout << "\n--- Back Substitution (U*x = y) ---\n";
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0;
                for (int j = i + 1; j < n; j++)
                    sum += U[i][j] * x[j];
                x[i] = (y[i] - sum) / U[i][i];
            }

            fout << "\nSolution Vector (x):\n";
            for (int i = 0; i < n; i++)
                fout << "x" << i + 1 << " = " << x[i] << "\n";

            // Verification: compute A*x and compare with b
            fout << "\n--- Verification (A*x = b) ---\n";
            vector<double> result(n);
            for (int i = 0; i < n; i++)
            {
                result[i] = 0;
                for (int j = 0; j < n; j++)
                    result[i] += A[i][j] * x[j];
                fout << "Row " << i + 1 << ": " << result[i] 
                     << " (expected: " << b[i] << ")\n";
            }
        }

        fout << "\n" << string(60, '=') << "\n\n";
    }

    fin.close();
    fout.close();
    cout << "All results written to output.txt\n";
    return 0;
}