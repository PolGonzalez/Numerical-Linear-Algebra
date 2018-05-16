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


// Decomposition A = LU.
void LU(int n, const Row& a, const Row& b, const Row& c, Row& beta, Row& alpha)  {
	alpha[0] = a[0];
	for (int i = 1; i < n; ++i) {
		beta[i] = b[i]/alpha[i - 1];
		alpha[i] = a[i] - beta[i]*c[i - 1];
	}
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
    
    Row q(n, 0);
    for (int i = 1; i < n; ++i) {
         double x;
         fileIn >> x;
         q[i] = x;
    }
    
    cout << "Dimensió del sistema: " << n << endl << endl;
    
	Row beta(n, 0);
	Row alpha(n, 0);
    LU(n, a, b, c, beta, alpha);
    
    
    // Hauríem de comprovar A - L*U.
    cout << " i  : Alpha  : Beta" << endl;
    for (int i = 0; i < n; ++i) {
        cout << i << ": " << alpha[i] << "  " << beta[i] << endl; 
    }
    
    Row x(n, 0);
    solution(n, beta, alpha, q, c, x);
    
    double epsilon = pow(10, -12);
    
    double norma = checker(a, b, c, n, beta, alpha, q, epsilon);
    cout << "Norma del màxim: " << norma << endl;
    
    norma = sub1(a, b, c, n, beta, alpha, q, epsilon);
    cout << "Norma_1: " << norma << endl;
    
    fileIn.close();
    
	
}
