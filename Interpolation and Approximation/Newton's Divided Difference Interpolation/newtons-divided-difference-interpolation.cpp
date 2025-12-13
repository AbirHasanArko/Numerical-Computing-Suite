#include <bits/stdc++.h>
using namespace std;

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
   Build Divided Difference Table
*/
vector<vector<double>> buildDividedDiffTable(const vector<double>& xs, const vector<double>& ys) {
    int n = (int)xs.size();
    vector<vector<double>> diff(n, vector<double>(n));
    
    // First column is y values
    for (int i=0; i<n; i++) diff[i][0] = ys[i];
    
    // Build divided difference table
    for (int j=1; j<n; j++){
        for (int i=0; i<n-j; i++){
            double denom = xs[i+j] - xs[i];
            if (fabs(denom) < 1e-15) {
                throw runtime_error("Duplicate x-values encountered");
            }
            diff[i][j] = (diff[i+1][j-1] - diff[i][j-1]) / denom;
        }
    }
    return diff;
}

/*
   Print Divided Difference Table
*/
void printDividedDiffTable(const vector<double>& xs, const vector<vector<double>>& diff, ostream& out) {
    int n = (int)xs.size();
    
    out << "\n====================================\n";
    out << "  DIVIDED DIFFERENCE TABLE\n";
    out << "====================================\n";
    out << setw(12) << "x" << setw(15) << "f(x)";
    for (int j=1; j<n; j++){
        out << setw(15) << ("f[x" + to_string(0) + "..." + "x" + to_string(j) + "]");
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
   Newton's Divided Difference Interpolation using pre-built table
*/
double newtonDividedDifferenceWithTable(const vector<double>& xs, const vector<vector<double>>& diff, double x) {
    int n = (int)xs.size();
    if (n == 1) return diff[0][0];
    
    double result = diff[0][0];
    double term = 1.0;
    
    for (int j=1; j<n; j++){
        term *= (x - xs[j-1]);
        result += diff[0][j] * term;
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
        double result = newtonDividedDifferenceWithTable(xs, diff, xInterpolate[i]);
        results[i] = result;
        
        bool isExtrap = (xInterpolate[i] < xs[0] || xInterpolate[i] > xs[n-1]);
        string extrapNote = isExtrap ? " (Extrapolation)" : "";
        
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
    printHeader("NEWTON'S DIVIDED DIFFERENCE INTERPOLATION", cout);
    
    string inputFile, outputFile;
    cout << "\nEnter input filename:  ";
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
    
    cout << "\nNumber of data points: " << n << "\n";
    
    // Open output file
    ofstream fout(outputFile);
    if (!fout) {
        cerr << "Error: Cannot create output file '" << outputFile << "'\n";
        return 1;
    }
    
    printHeader("NEWTON'S DIVIDED DIFFERENCE INTERPOLATION", fout);
    fout << "\nNumber of data points: " << n << "\n";
    
    // Print data points table
    printDataTable(xs, ys, cout);
    printDataTable(xs, ys, fout);
    
    // Build and print divided difference table
    try {
        vector<vector<double>> diffTable = buildDividedDiffTable(xs, ys);
        printDividedDiffTable(xs, diffTable, cout);
        printDividedDiffTable(xs, diffTable, fout);
        
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
            
            int nNew = n + 1;
            
            printHeader("  WITH ADDITIONAL DATA POINT", cout);
            printHeader("  WITH ADDITIONAL DATA POINT", fout);
            
            cout << "Additional point: x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            cout << "New number of data points: " << nNew << "\n";
            
            fout << "Additional point: x = " << fixed << setprecision(6) 
                 << xAdditional << ", y = " << yAdditional << "\n";
            fout << "New number of data points: " << nNew << "\n";
            
            // Print updated data points table
            printDataTable(xsNew, ysNew, cout);
            printDataTable(xsNew, ysNew, fout);
            
            // Build new divided difference table
            vector<vector<double>> diffTableNew = buildDividedDiffTable(xsNew, ysNew);
            printDividedDiffTable(xsNew, diffTableNew, cout);
            printDividedDiffTable(xsNew, diffTableNew, fout);
            
            // Re-interpolate and show differences
            printHeader("  UPDATED INTERPOLATION RESULTS", cout);
            printHeader("  UPDATED INTERPOLATION RESULTS", fout);
            
            for (int i=0; i<m; i++){
                double resultOld = interpolatedValues[i];
                double resultNew = newtonDividedDifferenceWithTable(xsNew, diffTableNew, xInterpolate[i]);
                double absError = fabs(resultNew - resultOld);
                double relError = (fabs(resultNew) > 1e-15) ? (absError / fabs(resultNew)) * 100.0 : 0.0;
                
                cout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                cout << "         Old y = " << setprecision(6) << resultOld << "\n";
                cout << "         New y = " << setprecision(6) << resultNew << "\n";
                cout << "         Absolute Difference: " << scientific << setprecision(6) << absError << "\n";
                cout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
                
                fout << "Point " << (i+1) << ": x = " << fixed << setprecision(6) << xInterpolate[i] << "\n";
                fout << "         Old y = " << setprecision(6) << resultOld << "\n";
                fout << "         New y = " << setprecision(6) << resultNew << "\n";
                fout << "         Absolute Difference: " << scientific << setprecision(6) << absError << "\n";
                fout << "         Relative Difference: " << fixed << setprecision(4) << relError << "%\n\n";
            }
            
            cout << "====================================\n";
            fout << "====================================\n";
        }
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        fout << "Error: " << e.what() << "\n";
        return 1;
    }
    
    fout.close();
    cout << "\nResults have been written to '" << outputFile << "'\n";
    cout << "Program executed successfully!\n";
    
    return 0;
}