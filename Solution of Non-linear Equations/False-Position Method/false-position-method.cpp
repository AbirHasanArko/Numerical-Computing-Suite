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

// Print polynomial in readable form
void printPolynomial(ofstream &out, const vector<f> &coef) {
    int n = coef.size() - 1;
    out << "f(x) = ";

    bool first = true;
    for (int i = 0; i <= n; i++) {
        f c = coef[i];
        int power = n - i;
        if (c == 0) continue;

        if (!first) {
            if (c > 0) out << " + ";
            else out << " - ";
        } else {
            if (c < 0) out << "-";
        }

        if (abs(c) != 1 || power == 0)
            out << abs(c);

        if (power > 0) {
            out << "x";
            if (power > 1)
                out << "^" << power;
        }

        first = false;
    }
    out << "\n\n";
}

// Print header
void printHeader(ofstream &out) {
    out << "============================================\n";
    out << "           False Position Method\n";
    out << "============================================\n\n";
}
void printHeader() {
    cout << "============================================\n";
    cout << "           False Position Method\n";
    cout << "============================================\n\n";
}

// Print roots table with trivial root highlight
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

int main() {
    string inputFile, outputFile;

    printHeader();
    // File names
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
    for (int i = 0; i <= degree; i++)
        in >> coef[i];

    printHeader(out);

    out << "Polynomial Degree : " << degree << "\n";
    out << "Coefficients      : ";
    for (auto c : coef) out << c << " ";
    out << "\n\n";

    printPolynomial(out, coef);

    // Parameters
    f searchRange = 5000.0, step = 0.5, tolerance = 1e-6;
    vector<f> roots, intervals;

    // Find intervals with sign changes
    for (f i = -searchRange; i <= searchRange; i += step) {
        f fx1 = fun(coef, i);
        f fx2 = fun(coef, i + step);

        // Exact root at grid
        if (abs(fx1) < tolerance) {
            bool duplicate = false;
            for (auto r : roots) if (abs(r - i) < tolerance) { duplicate = true; break; }
            if (!duplicate) roots.push_back(i);
        }
        // Sign change
        else if (fx1 * fx2 < 0.0) {
            intervals.push_back(i);
        }
    }

    // Apply False Position Method
    for (f start : intervals) {
        f xL = start, xR = start + step, x0;
        f fL = fun(coef, xL), fR = fun(coef, xR);

        // Make sure fL*fR < 0
        if (fL * fR > 0) continue;

        while (true) {
            // False Position formula
            x0 = (xL*fR - xR*fL) / (fR - fL);
            f f0 = fun(coef, x0);

            if (abs(f0) < tolerance || abs(xR - xL) < tolerance) {
                bool duplicate = false;
                for (auto r : roots) if (abs(r - x0) < tolerance) { duplicate = true; break; }
                if (!duplicate) roots.push_back(x0);
                break;
            }

            if (fL * f0 < 0) {
                xR = x0;
                fR = f0;
            } else {
                xL = x0;
                fL = f0;
            }
        }
    }

    sort(roots.begin(), roots.end());

    // Output roots table
    out << "Roots\n";
    out << "============================================\n";
    if (roots.empty()) {
        out << "No real roots found in the given range.\n\n";
    } else {
        printRootsTable(out, roots, coef, tolerance);
    }

    out << "============================================\n";
    out << "Computation Completed Successfully.\n";

    in.close();
    out.close();

    cout << "Computation completed. Results written to '" << outputFile << "'\n";

    return 0;
}
