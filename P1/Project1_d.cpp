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

    n = atoi(argv[1]);   // Dimension of the problem
    double h = 1./(1.+n);// Step length

    // Define arrays (vectors) to be used through Armadillo
    vec x = zeros<vec>(n+2);
    vec y_1 = zeros<vec>(n);
    vec u = zeros<vec>(n);
    vec diag = ones<vec>(n);
    // Calculate interval
    for (int i = 0; i <= n+1; i++){
        x[i] = i*h;
    }
    // Calculate RHS in the discrete Poisson's eq and analytical sol
    for (int i=1; i<=n; i++){
        y_1[i-1] = h*h*y(x[i]);
        u[i-1] = Sol(x[i]);

    }
    // Set up the tridiagonal matrix
    mat A = diagmat(2*diag);
    int j = 0;
    for (int i = 1; i<=n-2; i++){
        A(i,j) = -1;
        A(i,j+2) = -1;
        j += 1;
    }
    A(n-1,n-2) = -1;
    A(0,1) = -1;

    // Solve equations with LU-decomposition
    mat L, U;
    double start_t = clock();
    lu(L,U,A);

    vec b = solve(L,y_1);
    vec v = solve(U,b);

    // Calculate computation time
    double end_t = clock();
    double equation_time = (end_t - start_t)/(double)CLOCKS_PER_SEC;
    cout << "Equation time: " << equation_time << " seconds" << endl;

    //vec eps = log10(abs((u-v)/u));
    //cout << max(eps)<< endl;
    return 0;
}
