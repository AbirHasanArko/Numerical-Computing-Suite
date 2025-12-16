#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

int main() {
    int n;
    double a, b;
    cin >> n >> a >> b;

    vector<double> y(n+1);
    for(int i=0;i<=n;i++) cin >> y[i];

    if(n % 2 != 0){
        cout << "n must be even";
        return 0;
    }

    double h = (b-a)/n;
    double odd=0, even=0;

    for(int i=1;i<n;i++){
        if(i%2==0) even+=y[i];
        else odd+=y[i];
    }

    double result = (h/3)*(y[0]+y[n]+4*odd+2*even);
    cout << fixed << setprecision(6) << result;
    return 0;
}
