#include <bits/stdc++.h>
using namespace std;

static const double EPS = 1e-12;

static vector<vector<double>> getCofactor(const vector<vector<double>>& A, int p, int q) {
    int n = (int)A.size();
    vector<vector<double>> temp;
    temp.reserve(n - 1);

    for (int i = 0; i < n; i++) {
        if (i == p) continue;
        vector<double> row;
        row.reserve(n - 1);
        for (int j = 0; j < n; j++) {
            if (j == q) continue;
            row.push_back(A[i][j]);
        }
        temp.push_back(row);
    }
    return temp;
}

static double determinant(const vector<vector<double>>& A) {
    int n = (int)A.size();
    if (n == 1) return A[0][0];
    if (n == 2) return A[0][0] * A[1][1] - A[0][1] * A[1][0];

    double det = 0.0;
    int sign = 1;


    for (int f = 0; f < n; f++) {
        auto cof = getCofactor(A, 0, f);
        det += sign * A[0][f] * determinant(cof);
        sign = -sign;
    }
    return det;
}

static vector<vector<double>> cofactorMatrix(const vector<vector<double>>& A) {
    int n = (int)A.size();
    vector<vector<double>> C(n, vector<double>(n, 0.0));

    if (n == 1) {
        C[0][0] = 1.0;
        return C;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            auto cof = getCofactor(A, i, j);
            double minorDet = determinant(cof);
            double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            C[i][j] = sign * minorDet;
        }
    }
    return C;
}

static vector<vector<double>> transpose(const vector<vector<double>>& A) {
    int n = (int)A.size();
    vector<vector<double>> T(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            T[j][i] = A[i][j];
    return T;
}

static void printMatrix(const string& title, const vector<vector<double>>& A, int prec = 4) {
    cout << "\n" << title << "\n";
    cout << fixed << setprecision(prec);
    int n = (int)A.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(12) << A[i][j] << " ";
        }
        cout << "\n";
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cin >> n;
    if (!cin || n <= 0) {
        cout << "Invalid input.\n";
        return 0;
    }
    if (n > 8) {
        cout << "Warning: n is large for adjugate method (may be very slow). Recommended n <= 6.\n";
    }

    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            cin >> A[i][j];

    printMatrix("Input Matrix A:", A, 4);

    double detA = determinant(A);
    cout << "\nDeterminant det(A) = " << fixed << setprecision(6) << detA << "\n";

    if (fabs(detA) < EPS) {
        cout << "\nMatrix is singular (det(A) ≈ 0). Inverse does NOT exist.\n";
        return 0;
    }

    auto C = cofactorMatrix(A);
    auto adjA = transpose(C);

    printMatrix("Cofactor Matrix C:", C, 4);
    printMatrix("Adjugate Matrix adj(A) = C^T:", adjA, 4);

    vector<vector<double>> invA(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            invA[i][j] = adjA[i][j] / detA;

    printMatrix("Inverse Matrix A^{-1} = adj(A)/det(A):", invA, 6);

    return 0;
}
