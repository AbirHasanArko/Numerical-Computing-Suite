#include <bits/stdc++.h>
using namespace std;


static double f(int option, double x, double y) {
    switch (option) {
        case 1: return x + y;                 // y' = x + y
        case 2: return x - y;                 // y' = x - y
        case 3: return y - x*x + 1.0;         // y' = y - x^2 + 1
        case 4: return x * y;                 // y' = x*y
        default: return x + y;                // fallback
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int option;
    double x0, y0, xt, h;

    cin >> option;
    cin >> x0 >> y0;
    cin >> xt >> h;

    if (!cin || h <= 0) {
        cout << "Invalid input.\n";
        return 0;
    }
    if (xt < x0) {
        cout << "Target x must be >= x0.\n";
        return 0;
    }

    cout << fixed << setprecision(6);

    cout << "Runge-Kutta 4th Order (RK4)\n";
    cout << "Chosen function option = " << option << "\n";
    cout << "Initial condition: x0 = " << x0 << ", y0 = " << y0 << "\n";
    cout << "Target: x = " << xt << ", step size h = " << h << "\n\n";

    cout << left
         << setw(6)  << "Step"
         << setw(12) << "x"
         << setw(14) << "y"
         << setw(14) << "k1"
         << setw(14) << "k2"
         << setw(14) << "k3"
         << setw(14) << "k4"
         << "\n";

    cout << string(88, '-') << "\n";

    double x = x0, y = y0;
    int step = 0;


    int N = (int)round((xt - x0) / h);
    double x_end = x0 + N * h;

    if (fabs(x_end - xt) > 1e-9) {
        cout << "\nWarning: (xt - x0)/h is not an integer. "
             << "Program will stop at x = " << x_end << " (closest grid point).\n\n";
    }

    for (int i = 0; i < N; i++) {
        double k1 = h * f(option, x, y);
        double k2 = h * f(option, x + h/2.0, y + k1/2.0);
        double k3 = h * f(option, x + h/2.0, y + k2/2.0);
        double k4 = h * f(option, x + h,     y + k3);

        cout << left
             << setw(6)  << step
             << setw(12) << x
             << setw(14) << y
             << setw(14) << k1
             << setw(14) << k2
             << setw(14) << k3
             << setw(14) << k4
             << "\n";

        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        x = x + h;
        step++;
    }

    cout << "\nApproximate solution at x = " << x << " is y = " << y << "\n";
    return 0;
}
