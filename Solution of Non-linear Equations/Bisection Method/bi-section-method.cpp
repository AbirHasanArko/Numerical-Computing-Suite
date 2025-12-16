#include <bits/stdc++.h>
using namespace std;

#define f double

// Polynomial function f(x)
f fun(const vector<f> &coef, f x) {
    f result = 0.0;
    int n = coef.size() - 1;
    for (int i = 0; i <= n; i++)
        result += coef[i] * pow(x, n - i);
    return result;
}

// Print polynomial nicely
void printPolynomial(ofstream &out, const vector<f> &coef) {
    int n = coef.size() - 1;
    out << "f(x) = ";
    bool first = true;
    for (int i = 0; i <= n; i++) {
        f c = coef[i]; int power = n - i;
        if (c == 0) continue;
        if (!first) {
            if (c > 0) out << " + "; else out << " - ";
        } else { if (c < 0) out << "-"; }
        if (abs(c) != 1 || power == 0) out << abs(c);
        if (power > 0) { out << "x"; if (power > 1) out << "^" << power; }
        first = false;
    }
    out << "\n\n";
}

// Print header
void printHeader(ofstream &out) {
    out << "============================================\n";
    out << "              Bisection Method\n";
    out << "============================================\n\n";
}
void printHeader() {
    cout << "============================================\n";
    cout << "              Bisection Method\n";
    cout << "============================================\n\n";
}

// Print final roots table with trivial root highlight
void printRootsTable(ofstream &out, const vector<f> &roots, const vector<f> &coef, f tolerance = 1e-6) {
    out << setw(6) << "Index" 
        << setw(20) << "Root Value" << "\n";
    out << string(26, '-') << "\n";

    out << fixed << setprecision(6);
    for (size_t i = 0; i < roots.size(); i++) {
        f fx = fun(coef, roots[i]);
        out << setw(6) << (i + 1)
            << setw(20) << roots[i]
            << "\n";
    }
    out << "\n";
}

// Print iteration table per root
void printBisectionIteration(ofstream &out, const vector<tuple<int,f,f,f,f>> &iterData) {
    out << setw(10) << "Iteration"
        << setw(15) << "x_low"
        << setw(15) << "x_high"
        << setw(15) << "x_mid"
        << setw(20) << "f(x_mid)" << "\n";
    out << string(75, '-') << "\n";
    out << fixed << setprecision(6);
    for (auto &t : iterData) {
        int iter; f xL, xH, xM, fxM;
        tie(iter, xL, xH, xM, fxM) = t;
        out << setw(10) << iter
            << setw(15) << xL
            << setw(15) << xH
            << setw(15) << xM
            << setw(20) << fxM << "\n";
    }
    out << "\n";
}

int main() {
    string inputFile, outputFile;

    // File names
    printHeader();
    cout << "Enter input file name: ";
    cin >> inputFile;
    cout << "Enter output file name: ";
    cin >> outputFile;

    ifstream in(inputFile);
    ofstream out(outputFile);
    if (!in.is_open() || !out.is_open()) {
        cout << "Error: Cannot open input/output file!\n";
        return 1;
    }

    int degree;
    in >> degree;
    vector<f> coef(degree + 1);
    for (int i = 0; i <= degree; i++) in >> coef[i];

    printHeader(out);
    out << "Polynomial Degree : " << degree << "\n";
    out << "Coefficients      : ";
    for (auto c : coef) out << c << " ";
    out << "\n\n";
    printPolynomial(out, coef);

    // Bisection parameters
    f searchRange = 5000.0, step = 0.5, tolerance = 1e-6;
    vector<f> roots, intervals;

    // Scan for roots / intervals
    for (f i = -searchRange; i <= searchRange; i += step) {
        f fx1 = fun(coef, i);
        f fx2 = fun(coef, i + step);
        if (abs(fx1) < tolerance) {
            bool duplicate = false;
            for (auto r : roots) if (abs(r - i) < tolerance) { duplicate = true; break; }
            if (!duplicate) roots.push_back(i);
        } else if (fx1 * fx2 < 0.0) intervals.push_back(i);
    }

    vector<vector<tuple<int,f,f,f,f>>> allIterations; // iterations per root

    for (f start : intervals) {
        f xL = start, xH = start + step, xM;
        vector<tuple<int,f,f,f,f>> iterData;
        int iter = 1;
        while (true) {
            xM = (xL + xH) / 2.0;
            f fxM = fun(coef, xM);
            iterData.push_back({iter, xL, xH, xM, fxM});

            if (abs(fxM) < tolerance || abs(xH - xL) < tolerance) {
                bool duplicate = false;
                for (auto r : roots) if (abs(r - xM) < tolerance) { duplicate = true; break; }
                if (!duplicate) roots.push_back(xM);
                break;
            }

            f fxL = fun(coef, xL);
            if (fxL * fxM < 0) xH = xM;
            else xL = xM;

            iter++;
        }
        allIterations.push_back(iterData);
    }

    sort(roots.begin(), roots.end());

    // Print final roots
    out << "Roots\n";
    out << "============================================\n";
    if (roots.empty()) {
        out << "No real roots found in the given range.\n\n";
    } else {
        printRootsTable(out, roots, coef, tolerance);
    }

    // Print iteration tables
    for (size_t i = 0; i < allIterations.size(); i++) {
        out << "Iteration Table for Root " << (i + 1) << "\n";
        out << "----------------------------------------\n";
        printBisectionIteration(out, allIterations[i]);
    }

    out << "============================================\n";
    out << "Computation Completed Successfully.\n";

    in.close();
    out.close();
    cout << "Computation completed. Results written to '" << outputFile << "'\n";

    return 0;
}
