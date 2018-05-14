#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "stdlib.h"
using namespace std;

ifstream fileIn;
ofstream fileOut;

typedef vector<double> Row;
typedef vector<Row> Matrix;


double SupremeNorm(Matrix& A, int n){           // Norma del suprem matricial
    double max = 0;
    for (int i = 0; i < n; ++i){
        double sum = 0;
        for (int j = 0; j < n; ++j){
            sum += abs(A[i][j]);    
        }
        if (sum > max) max = sum;
    }
    return max;
}
double Norm1(Matrix& A, int n){                 // Norma sub 1 matricial
    double max = 0;
    for (int j = 0; j < n; ++j){
        double sum = 0;
        for (int i = 0; i < n; ++i){
            sum += abs(A[i][j]);    
        }
        if (sum > max) max = sum;
    }
    return max;
}
Matrix Product(Matrix& A, Matrix& B, int n){
    Matrix AB(n, Row(n, 0));

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            for (int k = 0 ; k < n; ++k){
                AB[i][j] += A[i][k]*B[k][j];
            }
        }
    } 
    return AB;
}
void PrintMatrix(Matrix& A, int n){
    for (int i = 0 ; i < n; ++i){
        for (int j = 0; j < n; ++j){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl << endl;
    return;
}


// Decomposition A = LU.
void LU(int n, const Row& a, const Row& b, const Row& c, double d, double e, Row& beta, Row& alpha, Row& epsilon, Row& delta)  {
    // alpha i beta
	alpha[0] = a[0];
	for (int i = 1; i < n - 1; ++i) {
		beta[i] = b[i]/alpha[i - 1];
		alpha[i] = a[i] - beta[i]*c[i - 1];
	}
    
    // epsilon
    epsilon[0] = e;
    for (int i = 1; i < n - 1; ++i) {
        epsilon[i] = -beta[i]*epsilon[i - 1];
    }
    epsilon[n - 2] = c[n - 2] - beta[n - 2]*epsilon[n - 3];
    
    // delta
    delta[0] = d/alpha[0];
    for (int i = 1; i < n - 2; ++i) {
        delta[i] = (-delta[i - 1]*c[i - 1])/alpha[i];
    }
    delta[n - 2] = (b[n - 1] - delta[n - 3]*c[n - 3])/alpha[n - 2];
    
    double sum = 0;
    for (int i = 0; i < n - 1; ++i) {
        sum += epsilon[i]*delta[i];
    }
    alpha[n - 1] = a[n - 1] - sum;
}

// Solves Ly = b and Ux = y.
// Solves Ly = b and Ux = y.
void solution(int n, const Row& beta, const Row& alpha, const Row& p, const Row& c, Row& x) {
	Row y(n, 0);
	y[0] = p[0];
	for (int i = 1; i < n; ++i) {
		y[i] = p[i] - beta[i]*y[i - 1];
	}
	
	
	x[n - 1] = y[n - 1]/alpha[n - 1];
	for (int i = n - 2; i >= 0; --i) {
		x[i] = (y[i] - x[i + 1]*c[i])/alpha[i];
	}
}

// Returns true if the solution obtained has an error smaller than epsilon.
double checker(const Row& a, const Row& b, const Row& c, int n, const Row& beta, const Row& alpha, const Row& p, double epsilon) {
	Row R(n, 0);
	R[0] = a[0]*p[0] + c[0]*p[1] - p[0];
	for (int i = 1; i < n; ++i) {
		R[i] = b[i]*p[i - 1] + a[i]*p[i] + c[i]*p[i + 1] -p[i];
	}
	R[n - 1] = b[n - 1]*p[n - 2] + a[n - 1]*p[n - 1] - p[n - 1];
	
	double max = abs(R[0]);
	for (int i = 1; i < n; ++i) {
		if (abs(R[i]) > max) max = abs(R[i]);
	}
	
	return max;
}

double sub1(const Row& a, const Row& b, const Row& c, int n, const Row& beta, const Row& alpha, const Row& p, double epsilon) {
	Row R(n, 0);
	R[0] = a[0]*p[0] + c[0]*p[1] - p[0];
	for (int i = 1; i < n; ++i) {
		R[i] = b[i]*p[i - 1] + a[i]*p[i] + c[i]*p[i + 1] -p[i];
	}
	R[n - 1] = b[n - 1]*p[n - 2] + a[n - 1]*p[n - 1] - p[n - 1];
	
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		sum += abs(R[i]);
	}
	
	return sum;
}

void matrius(Row& c, Row& beta, Row& alpha, Row& epsilon, Row& delta, Matrix& L, Matrix& U) {
    int n = (int)beta.size();
    // L
    for (int i = 0; i < n; ++i) L[i][i] = 1;
    for (int i = 1; i < n - 1; ++i) L[i][i - 1] = beta[i];
    for (int i = 0; i < n - 1; ++i) L[n - 1][i] = delta[i];
    cout << "Matriu L: " << endl;
    PrintMatrix(L, n);
    
    // U
    for (int i = 0; i < n; ++i) U[i][i] = alpha[i];
    for (int i = 0; i < n - 2; ++i) U[i][i + 1] = c[i];
    for (int i = 0; i < n - 1; ++i) U[i][n - 1] = epsilon[i];
    cout << "Matriu U: " << endl;
    PrintMatrix(U, n);
}


int main(int argc, char *argv[]) {
    
	if (argc >= 2) fileIn.open(argv[1]); // obrim el fitxer per llegir les dades.
	else {cerr << "ERROR: no hi ha fitxer d'entrada" << endl; return 0;}
    
    if (!fileIn.is_open()) {
        cerr << "Error: no s'ha pogut obrir el fitxer " << argv[1]
            << " per llegir les dades" << endl;
        return 0;
    }
    
    int n;
    fileIn >> n;
    
    Row a(n, 0);
    for (int i = 0; i < n; ++i) {
        double x;
        fileIn >> x;
        a[i] = x;
    }
    
    Row b(n, 0);
    for (int i = 1; i < n; ++i) {
         double x;
         fileIn >> x;
         b[i] = x;
    }
    
    Row c(n, 0);
    for (int i = 0; i < n - 1; ++i) {
        double x;
        fileIn >> x;
        c[i] = x;
    }
    
    /*double d;
    cin >> d;
    
    double e;
    cin >> e;*/
    double d = 1, e = -1;
    
    Row q(n, 0);
    for (int i = 1; i < n; ++i) {
         double x;
         fileIn >> x;
         q[i] = x;
    }
    
    fileIn.close();
    
    cout << "Dimensió del sistema: " << n << endl << endl;
    
	Row beta(n, 0);
	Row alpha(n, 0);
    Row epsilon(n, 0);
    Row delta (n, 0);
    LU(n, a, b, c, d, e, beta, alpha, epsilon, delta);
    
    // comprovem la norma A - L*U construint les matrius L i U
    Matrix L(n, Row(n, 0));
    Matrix U(n, Row(n, 0));
    matrius(c, beta, alpha, epsilon, delta, L, U);
    
    
    Matrix res = Product(L, U, n);
    res[n - 1][0] -= d;
    res[0][n - 1] -= e;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i][i] -= a[i];
            if (j != n - 1) { 
                res[i][i + 1] -= c[i];
                res[i + 1][j] -= b[i];
            }
        }
    }
    
    cout << "norma sub 1 matricial: " << Norm1(res, n) << endl;
    cout << "norma sub infinit matricial: " << SupremeNorm(res, n) << endl;
    
    
    cout << " i  : Alpha  : Beta : Delta : Epsilon" << endl;
    for (int i = 0; i < n; ++i) {
        cout << i << ": " << alpha[i] << "  " << 
        beta[i] << " " << delta[i] << " " << epsilon[i] << endl; 
    }
	
}
