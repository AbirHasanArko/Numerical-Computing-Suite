#include <bits/stdc++.h>
using namespace std;

/*
   Utility:  factorial up to n (double)
*/
double factorial(int n){
    double f = 1.0;
    for (int i=2; i<=n; i++) f *= i;
    return f;
}

/*
   Print Data Points Table
*/
void printDataTable(const vector<double>& xs, const vector<double>& ys, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  DATA POINTS TABLE\n";
    out << "====================================\n";
    
    // Print x values
    out << setw(10) << "x";
    for (int i=0; i<n; i++){
        out << setw(15) << fixed << setprecision(4) << xs[i];
    }
    out << "\n";
    
    // Print separator
    out << string(10 + 15*n, '-') << "\n";
    
    // Print y values
    out << setw(10) << "y";
    for (int i=0; i<n; i++){
        out << setw(15) << fixed << setprecision(6) << ys[i];
    }
    out << "\n";
    out << "====================================\n";
}

/*
   Build Forward Difference Table
*/
vector<vector<double>> buildForwardDiffTable(const vector<double>& xs, const vector<double>& ys) {
    int n = (int)xs.size();
    vector<vector<double>> diff(n, vector<double>(n));
    
    for (int i=0; i<n; i++) diff[i][0] = ys[i];
    for (int j=1; j<n; j++){
        for (int i=0; i<n-j; i++){
            diff[i][j] = diff[i+1][j-1] - diff[i][j-1];
        }
    }
    return diff;
}

/*
   Print Forward Difference Table
*/
void printForwardDiffTable(const vector<double>& xs, const vector<vector<double>>& diff, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  FORWARD DIFFERENCE TABLE\n";
    out << "====================================\n";
    out << setw(12) << "x" << setw(15) << "y";
    for (int j=1; j<n; j++){
        out << setw(15) << ("Delta^" + to_string(j) + "y");
    }
    out << "\n" << string(12 + 15*n, '-') << "\n";
    
    for (int i=0; i<n; i++){
        out << fixed << setprecision(4) << setw(12) << xs[i];
        for (int j=0; j<n-i; j++){
            out << setprecision(6) << setw(15) << diff[i][j];
        }
        out << "\n";
    }
    out << "====================================\n";
}

/*
   Newton's Forward Interpolation using pre-built table
*/
double newtonForwardWithTable(const vector<double>& xs, const vector<vector<double>>& diff, double x) {
    int n = (int)xs.size();
    if (n == 1) return diff[0][0];
    
    double h = xs[1] - xs[0];
    double u = (x - xs[0]) / h;
    double result = diff[0][0];
    double u_prod = 1.0;
    
    for (int k=1; k<n; k++){
        u_prod *= (u - (k-1));
        result += (u_prod / factorial(k)) * diff[0][k];
    }
    return result;
}

/*
   Process and output interpolation results
*/
void processInterpolation(const vector<double>& xs, const vector<vector<double>>& diff,
                         const vector<double>& xInterpolate, vector<double>& results,
                         ostream& cout_stream, ostream& fout) {
    int m = xInterpolate.size();
    int n = xs.size();
    
    for (int i=0; i<m; i++){
        double result = newtonForwardWithTable(xs, diff, xInterpolate[i]);
        results[i] = result;
        
        bool isExtrap = (xInterpolate[i] < xs[0] || xInterpolate[i] > xs[n-1]);
        string extrapNote = isExtrap ?  " (Extrapolation)" : "";
        
        string output = "Point " + to_string(i+1) + ": x = ";
        cout_stream << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        cout_stream << "         y = " << setprecision(6) << result << "\n";
        
        fout << output << fixed << setprecision(6) << xInterpolate[i] << extrapNote << "\n";
        fout << "         y = " << setprecision(6) << result << "\n";
    }
}

/*
   Print section header
*/
void printHeader(const string& title, ostream& out) {
    out << "\n====================================\n";
    out << title << "\n";
    out << "====================================\n";
}

/*
   Main with File I/O
*/
int main() {
    printHeader("NEWTON'S FORWARD INTERPOLATION", cout);
    
    string inputFile, outputFile;
    cout << "\nEnter input filename: ";
    cin >> inputFile;
    cout << "Enter output filename: ";
    cin >> outputFile;
    
    // Open input file
    ifstream fin(inputFile);
    if (!fin) {
        cerr << "Error: Cannot open input file '" << inputFile << "'\n";
        return 1;
    }
    
    // Read data points
    int n;
    fin >> n;
    if (n <= 0) {
        cerr << "Error: Invalid number of data points\n";
        return 1;
    }
    
    vector<double> xs(n), ys(n);
    for (int i=0; i<n; i++){
        fin >> xs[i] >> ys[i];
    }
    
    // Read interpolation points
    int m;
    fin >> m;
    if (m <= 0) {
        cerr << "Error: Invalid number of interpolation points\n";
        return 1;
    }
    
    vector<double> xInterpolate(m);
    for (int i=0; i<m; i++){
        fin >> xInterpolate[i];
    }
    
    // Check for additional data point
    double xAdditional, yAdditional;
    bool hasAdditional = false;
    if (fin >> xAdditional >> yAdditional) {
        hasAdditional = true;
    }
    fin.close();
    
    // Verify equal spacing
    if (n > 1) {
        double h = xs[1] - xs[0];
        for (int i=2; i<n; i++){
            if (fabs((xs[i] - xs[i-1]) - h) > 1e-9) {
                cerr << "Error:  Data points are not equally spaced!\n";
                return 1;
            }
        }
        cout << "\nNumber of data points: " << n << "\n";
        cout << "Step size (h): " << fixed << setprecision(6) << h << "\n";
    }
    
    // Open output file
    ofstream fout(outputFile);
    if (!fout) {
        cerr << "Error:  Cannot create output file '" << outputFile << "'\n";
        return 1;
    }
    
    printHeader("NEWTON'S FORWARD INTERPOLATION", fout);
    fout << "\nNumber of data points: " << n << "\n";
    if (n > 1) {
        fout << "Step size (h): " << fixed << setprecision(6) << (xs[1] - xs[0]) << "\n";
    }
    
    // Print data points table
    printDataTable(xs, ys, cout);
    printDataTable(xs, ys, fout);
    
    // Build and print forward difference table
    vector<vector<double>> diffTable = buildForwardDiffTable(xs, ys);
    printForwardDiffTable(xs, diffTable, cout);
    printForwardDiffTable(xs, diffTable, fout);
    
    // Perform interpolation
    printHeader("  INTERPOLATION RESULTS", cout);
    printHeader("  INTERPOLATION RESULTS", fout);
    
    vector<double> interpolatedValues(m);
    processInterpolation(xs, diffTable, xInterpolate, interpolatedValues, cout, fout);
    
    cout << "====================================\n";
    fout << "====================================\n";
    
    // Process additional point if provided
    if (hasAdditional) {
        vector<double> xsNew = xs;
        vector<double> ysNew = ys;
        
        // Insert in sorted order
        int insertPos = n;
        for (int i=0; i<n; i++){
            if (xAdditional < xs[i]) {
                insertPos = i;
                break;
            }
        }
        xsNew.insert(xsNew.begin() + insertPos, xAdditional);
        ysNew.insert(ysNew. begin() + insertPos, yAdditional);
        
        // Verify equal spacing
        int nNew = n + 1;
        double hNew = xsNew[1] - xsNew[0];
        bool validSpacing = true;
        for (int i=2; i<nNew; i++){
            if (fabs((xsNew[i] - xsNew[i-1]) - hNew) > 1e-9) {
                validSpacing = false;
                break;
            }
        }
        
        if (validSpacing) {
            printHeader("  WITH ADDITIONAL DATA POINT", cout);
            printHeader("  WITH ADDITIONAL DATA POINT", fout);
            
            cout << "Additional point:  x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            cout << "New number of data points: " << nNew << "\n";
            cout << "New step size (h): " << setprecision(6) << hNew << "\n";
            
            fout << "Additional point: x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            fout << "New number of data points: " << nNew << "\n";
            fout << "New step size (h): " << setprecision(6) << hNew << "\n";
            
            // Print updated data points table
            printDataTable(xsNew, ysNew, cout);
            printDataTable(xsNew, ysNew, fout);
            
            // Build new difference table
            vector<vector<double>> diffTableNew = buildForwardDiffTable(xsNew, ysNew);
            printForwardDiffTable(xsNew, diffTableNew, cout);
            printForwardDiffTable(xsNew, diffTableNew, fout);
            
            // Re-interpolate and show differences
            printHeader("  UPDATED INTERPOLATION RESULTS", cout);
            printHeader("  UPDATED INTERPOLATION RESULTS", fout);
            
            for (int i=0; i<m; i++){
                double resultOld = interpolatedValues[i];
                double resultNew = newtonForwardWithTable(xsNew, diffTableNew, xInterpolate[i]);
                double absError = fabs(resultNew - resultOld);
                double relError = (fabs(resultNew) > 1e-15) ? (absError / fabs(resultNew)) * 100.0 : 0.0;
                
                cout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                cout << "         Old y = " << setprecision(6) << resultOld << "\n";
                cout << "         New y = " << setprecision(6) << resultNew << "\n";
                cout << "         Absolute Difference:   " << scientific << setprecision(6) << absError << "\n";
                cout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
                
                fout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                fout << "         Old y = " << setprecision(6) << resultOld << "\n";
                fout << "         New y = " << setprecision(6) << resultNew << "\n";
                fout << "         Absolute Difference:  " << scientific << setprecision(6) << absError << "\n";
                fout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
            }
            
            cout << "====================================\n";
            fout << "====================================\n";
        } else {
            printHeader("WARNING: Additional point breaks equal spacing!", cout);
            cout << "Cannot use Newton's Forward Interpolation.\n";
            cout << "====================================\n";
            
            printHeader("WARNING: Additional point breaks equal spacing!", fout);
            fout << "Cannot use Newton's Forward Interpolation.\n";
            fout << "====================================\n";
        }
    }
    
    fout.close();
    cout << "\nResults have been written to '" << outputFile << "'\n";
    cout << "Program executed successfully!\n";
    
    return 0;
}