#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include<ctime>

using namespace std;
using namespace arma;
ofstream fs;
// Analytical solution for u(x)
double Sol(double x){
    return 1. + (exp(-10.) -1.)*x - exp(-10.*x);
}
// RHS in the poisson equation
double y(double x){
    return 100.*exp(-10.*x);
}
// Code for general defined values in array a,b,c given {a}={c} = -1 and {b} = b = 2
int main(int argc, char* argv[]){
    int n;
    char *ofilename;
    if(argc < 2){
        cout << endl << "Argument not given to" << argv[0] << endl <<
                "Set an integer 'n'" << endl << "To save values include a file name: filename.txt (optional)" << endl;
        exit(1);
    }

    n = atoi(argv[1]);

    int b = 2;  // b values in the diagonal, b_2 are the alterated, see below.


    double h = 1./(n + 1.);        // Step width in x
    double *v = new double[n+2];   // To contain numerical solution
    double *u = new double[n+2];   // To contain analytical solution
    double *x = new double[n+2];   // To contain interval x in [0,1] (but use interval (0,1) )
    double *y_1 = new double[n+1]; // To contain the RHS in descritized algorithm
    double *b_2 = new double[n+1]; // To contain new b_2 values used in forward substitution
    double *y_2 = new double[n+1]; // To contain new y_2 values used in --"--
    // Set values and initial conditions, bc's
    u[0] = 0.;
    v[0] = 0.;
    for(int i=0; i<=n+1; i++){
        x[i] = i*h;
    }
    for(int i = 1; i<=n; i++){
        y_1[i] = h*h*y(x[i]);
        u[i] = Sol(x[i]);
    }
    b_2[1] = b;
    y_2[1] = y_1[1];
    // Forward substitution, calculate values in order to find num sol v
    double start_t = clock();
    for(int i=2; i <= n; i++){
        b_2[i] = (i + 1.)/i;
        y_2[i] = y_1[i] + y_2[i-1]/b_2[i-1]; // can check testfile by changing + to - and run test
    }

    v[n] = y_2[n]/b_2[n]; // Set last value v(n)
    for(int i=n-1; i>= 1; i--){
        v[i] = (y_2[i] + v[i+1])/b_2[i];  // calculating numerical solution
    }
    // Calculate equation time
    double end_t = clock();
    double equation_time = (end_t - start_t)/(double)CLOCKS_PER_SEC;
    cout << "Equation time: " << equation_time << " seconds" << endl;

    // Write data to file (or if filename = Test.txt, run test)
    if (argc == 3){
        ofilename = argv[2];
        string test = "Test.txt";
        bool equal = (ofilename == test);

        if (equal == 1){  // If user is running test
            if (n != 1000){
                cout << "Set n = 1000" << endl;
                exit(2);
            }
            else{
                fs.open(ofilename);
                cout << "Creating temporary test file.." << endl;
                fs << setiosflags(ios::showpoint | ios::uppercase);
                std::setprecision(8);
                for (int i = 1; i <= n; i++){
                    fs << x[i] << '\t' << u[i] << '\t' << v[i] << endl;
                }
                fs.close();
                cout << "Test file created.." << endl;
                // Call python, run the test there (see: relerror.py lines 31-53)
                string pyscript = " relerror.py Test.txt";
                string command = "python";
                command += pyscript;
                system(command.c_str());
            }

        }
        else { // Write values to file
            string input;
            cout << "Want to save values in file: "<< argv[2] << "? y/n" << endl;
            cin >> input;
            if (input == "yes" || input == "y"){
                fs.open(ofilename);
                fs << setiosflags(ios::showpoint | ios::uppercase);
                cout << endl << "writing values in columns: x, u, v" << endl << "Column size: "
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


    }
    else {
        cout << endl << "Values were not stored, include a filename in order to store values: filename.txt" << endl;
    }
    // Free from memory
    delete [] x;
    delete [] v;
    delete [] u;
    delete [] y_1;
    delete [] y_2;
    delete [] b_2;

    return 0;
}
