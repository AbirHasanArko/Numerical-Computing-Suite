#include <bits/stdc++.h>
using namespace std;

double determinant(const vector<vector<double>>& A);
vector<vector<double>> adjoint(const vector<vector<double>>& A);

/* ---------------------------
   Compute inverse by 1/det(A) * adj(A)
----------------------------*/
bool inverseByAdjoint(const vector<vector<double>>& A,
                      vector<vector<double>>& inv) {
    int n = A.size();
    for (auto &r : A)
        if (r.size() != n) return false; // must be square

    double detA = determinant(A);
    if (fabs(detA) < 1e-12) return false; // non-invertible

    vector<vector<double>> adjA = adjoint(A);

    inv.assign(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adjA[i][j] / detA;

    return true;
}

/* ---------------------------
   Determinant (cofactor expansion)
----------------------------*/
double determinant(const vector<vector<double>>& A) {
    int n = A.size();
    if (n == 1) return A[0][0];
    if (n == 2)
        return A[0][0]*A[1][1] - A[0][1]*A[1][0];

    double det = 0;
    for (int col = 0; col < n; col++) {
        vector<vector<double>> sub(n-1, vector<double>(n-1));
        for (int i = 1; i < n; i++) {
            int c2 = 0;
            for (int j = 0; j < n; j++) {
                if (j == col) continue;
                sub[i-1][c2++] = A[i][j];
            }
        }
        det += ((col % 2 == 0) ? 1 : -1) * A[0][col] * determinant(sub);
    }
    return det;
}

/* ---------------------------
   Adjoint matrix
----------------------------*/
vector<vector<double>> adjoint(const vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> adj(n, vector<double>(n));

    if (n == 1) {
        adj[0][0] = 1;
        return adj;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vector<vector<double>> sub(n-1, vector<double>(n-1));
            int r2 = 0;
            for (int r = 0; r < n; r++) {
                if (r == i) continue;
                int c2 = 0;
                for (int c = 0; c < n; c++) {
                    if (c == j) continue;
                    sub[r2][c2++] = A[r][c];
                }
                r2++;
            }

            double cofactor = ((i + j) % 2 == 0 ? 1 : -1)
                              * determinant(sub);
            adj[j][i] = cofactor; // transpose
        }
    }
    return adj;
}

/* ---------------------------
   Main: File I/O
----------------------------*/
int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin.is_open()) {
        cerr << "Failed to open input.txt\n";
        return 1;
    }

    int n;
    fin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fin >> A[i][j];

    vector<vector<double>> inv;
    fout << fixed << setprecision(6);

    if (!inverseByAdjoint(A, inv)) {
        fout << "Matrix is singular.\n";
        return 0;
    }

    fout << "Inverse using adj(A)/det(A):\n";
    for (auto &row : inv) {
        for (double x : row)
            fout << setw(12) << x;
        fout << "\n";
    }

    fin.close();
    fout.close();
    return 0;
}
