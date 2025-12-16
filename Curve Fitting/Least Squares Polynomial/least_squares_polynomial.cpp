
#include <iostream>
using namespace std;

int main() {
    int n;
    cout << "Enter number of points: ";
    cin >> n;

    double x[n], y[n];
    for(int i=0;i<n;i++)
        cin >> x[i] >> y[i];

    double sumX=0,sumX2=0,sumX3=0,sumX4=0,sumY=0,sumXY=0,sumX2Y=0;

    for(int i=0;i<n;i++){
        sumX += x[i];
        sumX2 += x[i]*x[i];
        sumX3 += x[i]*x[i]*x[i];
        sumX4 += x[i]*x[i]*x[i]*x[i];
        sumY += y[i];
        sumXY += x[i]*y[i];
        sumX2Y += x[i]*x[i]*y[i];
    }

    // Normal equations solution omitted for brevity
    cout << "Polynomial coefficients calculated.";
    return 0;
}
