#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <ctime>

using namespace std;
using namespace arma;
ofstream fs;
// Analytical solution u(x)
double Sol(double x){
    return 1. + (exp(-10.) -1.)*x - exp(-10.*x);
}
// RHS in the poisson's equation
double y(double x){
    return 100.*exp(-10.*x);
}
// Code for general defined values in array a,b,c not necessary equal
int main(int argc, char* argv[]){
    int n;
    char *ofilename;
    if(argc < 2){
        cout << endl << "Argument not given to" << argv[0] << endl <<
                "Set an integer 'n'" << endl << "To save values include a file name: filename.txt (optional)" << endl;
        exit(1);
    }

    n = atoi(argv[1]);
    // Define arrays = diagonal('s') in the tridiagonal matrix
    int *a = new int[n+1];
    int *b = new int[n+1];
    int *c = new int[n+1];

    double h = 1./(n + 1.);      // Step width in x
    double *v = new double[n+2]; // To contain numerical solution
    double *u = new double[n+2]; // To contain analytical solution
    double *x = new double[n+2]; // To contain interval x in (0,1)
    double *y_1 = new double[n+1]; // To contain the RHS in descritized algorithm
    double *b_2 = new double[n+1]; // To contain new b_2 values used in forward substitution
    double *y_2 = new double[n+1]; // To contain new y_2 values used in --"--
    // Set values and initial conditions and interval
    u[0] = 0.;
    v[0] = 0.;
    for(int i=0; i<=n+1; i++){
        x[i] = i*h;
    }
    for(int i = 1; i<=n; i++){
        y_1[i] = h*h*y(x[i]);
        u[i] = Sol(x[i]);
        a[i] = -1.;
        c[i] = -1.;
        b[i] = 2.;
    }
    a[1] = 0.;
    c[n] = 0.;
    b_2[1] = b[1];
    y_2[1] = y_1[1];
    // Forward substitution, calculate values in order to find num sol v
    double start_t = clock();
    for(int i=2; i <= n; i++){
        b_2[i] = b[i] - a[i]*c[i-1]/b_2[i-1];
        y_2[i] = y_1[i] - a[i]*y_2[i-1]/b_2[i-1];
    }

    v[n] = y_2[n]/b_2[n]; // Set last value v(n) and iterate
    for(int i=n-1; i>= 1; i--){
        v[i] = (y_2[i] - c[i]*v[i+1])/b_2[i];
    }
    double end_t = clock();
    double equation_time = (end_t - start_t)/(double)CLOCKS_PER_SEC;
    cout << endl << "Equation time: " << equation_time << " seconds" << endl;
    //cout << v[-1] << " " << u[-1];

    // Save values in text file
    if (argc == 3){
        string input;
        ofilename = argv[2];
        cout << "Want to save values in file: "<< argv[2] << "? y/n" << endl; // This is not necessary, but easier to check code
        cin >> input;
        if (input == "yes" || input == "y"){
            fs.open(ofilename);
            fs << setiosflags(ios::showpoint | ios::uppercase);
            cout << endl << "saving values in columns: x, u, v" << endl << "Column size: "
                 << n << endl;
            std::setprecision(8);
            for (int i = 1; i <= n; i++){
                fs << x[i] << '\t' << u[i] << '\t' << v[i] << endl;
            }
            fs.close();
            cout << endl << "Values stored." << endl;
        }
        else{
            cout << endl << "Values not stored." << endl;
            }

    }
    else {
        cout << endl << "Values were not stored, include a filename in order to store values: filename.txt" << endl;
    }
    // Free from memory
    delete [] x;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] v;
    delete [] u;
    delete [] y_1;
    delete [] y_2;
    delete [] b_2;

    return 0;
}
