#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "stdlib.h"
using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix;

double norminf(const Row& v, int n) {						// Norma sub infinit vectorial
	double maxim = -1;
	for (int i = 0; i < n; ++i) {
		maxim = max(maxim, abs(v[i]));
	}
	return maxim;
}
double norma1(const Row& v, int n) {						// Norma sub 1 vectorial
    double sum = 0; 
    for (int i = 0; i < n; ++i) {
        sum += abs(v[i]);
    }
    return sum;
}
double SupremeNorm(const Matrix& A, int n){          		// Norma del suprem matricial
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
double Norm1(const Matrix& A, int n){               	 	// Norma sub 1 matricial
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
Matrix Product(const Matrix& A, Matrix& B, int n){			// Producte de dues matriu
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
void PrintMatrix(const Matrix& A, int n){					// Escriu una matriu
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
void solution(int n, const Row& c, const Row& beta, const Row& alpha, const Row& epsilon, const Row& delta, const Row& b, Row& x) {
	// Ly = b
	Row y(n, 0);	
	y[0] = b[0];
	for (int i = 1; i < n - 1; ++i) {
		y[i] = b[i] - beta[i]*y[i - 1];
	}
	double sum = 0;
	for (int i = 0; i < n - 1; ++i) sum += delta[i]*y[i];
	y[n - 1] = b[n - 1] - sum;
	
	// Ux = y
	x[n - 1] = y[n - 1]/alpha[n - 1];
	x[n - 2] = (y[n - 2] - epsilon[n - 2]*x[n - 1])/alpha[n - 2];
	for (int i = n - 3; i >= 0; --i) {
		x[i] = (y[i] - epsilon[i]*x[n - 1] - c[i]*x[i + 1])/alpha[i];
	}
}

// Construeix les matrius L i U
void matrius(const Row& c, const Row& beta, const Row& alpha, const Row& epsilon, const Row& delta, Matrix& L, Matrix& U) {
    int n = (int)beta.size();
    // L
    for (int i = 0; i < n; ++i) L[i][i] = 1;
    for (int i = 1; i < n - 1; ++i) L[i][i - 1] = beta[i];
    for (int i = 0; i < n - 1; ++i) L[n - 1][i] = delta[i];
    // U
    for (int i = 0; i < n; ++i) U[i][i] = alpha[i];
    for (int i = 0; i < n - 2; ++i) U[i][i + 1] = c[i];
    for (int i = 0; i < n - 1; ++i) U[i][n - 1] = epsilon[i];
}


int main(int argc, char *argv[]) {
	
	ifstream fileIn;
    ofstream fileOut;
    
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
    
    cout << "System dimension: " << n << endl << endl;
    
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
    	res[i][i] -= a[i];
    	if (i != n - 1) res[i][i + 1] -= c[i];
    	if (i != n - 1) res[i + 1][i] -= b[i + 1];
	}
	
    cout << "Descomposicio A = LU: OK" << endl;
    cout << "Norma del residu L*U - A." << endl;
    cout << "   norma sub 1 matricial: " << Norm1(res, n) << endl;
    cout << "   norma sub infinit matricial: " << SupremeNorm(res, n) << endl << endl;
    cout << "Coeficients de la descomposicio: " << endl;
    cout << " i  : Alpha  : Beta : Delta : Epsilon" << endl;
    for (int i = 0; i < n; ++i) {
        cout << i << ": " << alpha[i] << "  " << 
        beta[i] << " " << delta[i] << " " << epsilon[i] << endl; 
    }
    
    // Resolució del sistema
    Row x(n, 0);
    solution(n, c, beta, alpha, epsilon, delta, q, x);
    
    Row r(n, 0);
    r[0] = a[0]*x[0] + c[0]*x[1] + e*x[n - 1];
    for (int i = 1; i < n - 1; ++i) r[i] = b[i]*x[i - 1] + a[i]*x[i] + c[i]*x[i + 1];
    r[n - 1] = d*x[0] + b[n - 1]*x[n - 2] + a[n - 1]*x[n - 1];
    
    for (int i = 0; i < n; ++i) {
    	r[i] -= q[i];
	}
	cout << endl;
	cout << "Resolucio del sistema Ax = q: OK" << endl;
	cout << "Norma del residu A*x - q" << endl;
	cout << "    norma sub 1 vectorial: " << norma1(r, n) << endl;
	cout << "    norma sub infinit vectorial: " << norminf(r, n) << endl;
	cout << "Vector solucio: " << endl;
	cout << "i  :  x" << endl;
	for (int i = 0; i < n; ++i) {
		cout << i << ": " << x[i] << endl;
	}
}
