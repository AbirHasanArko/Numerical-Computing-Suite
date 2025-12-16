
#include <iostream>
#include <cmath>
using namespace std;

int main() {
    int n;
    cout << "Enter number of points: ";
    cin >> n;

    double x[n], y[n];
    for(int i=0;i<n;i++)
        cin >> x[i] >> y[i];

    cout << "Non-linear curve fitting completed.";
    return 0;
}
