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

    if(n % 3 != 0){
        cout << "n must be multiple of 3";
        return 0;
    }

    double h = (b-a)/n;
    double sum3=0, sum2=0;

    for(int i=1;i<n;i++){
        if(i%3==0) sum2+=y[i];
        else sum3+=y[i];
    }

    double result = (3*h/8)*(y[0]+y[n]+3*sum3+2*sum2);
    cout << fixed << setprecision(6) << result;
    return 0;
}
