#include <iostream>
#include <cmath>
using namespace std;


// Diferencia A-B, amb A i B matrius nxn
double **Difference(double** A, double** B, int n){
    double **D;
    D = (double**) calloc(n, sizeof(double*));
    for (int i = 0; i < n; ++i) D[i] = (double *) calloc(n, sizeof(double));

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
                D[i][j] = A[i][j] - B[i][j];
            }
    }
    return D;
}

// Multiplica dues matrius nxn
double **Product(double** A, double** B, int n){
    double **AB;
    AB = (double**) calloc(n, sizeof(double*));
    for (int i = 0; i < n; ++i) AB[i] = (double *) calloc(n, sizeof(double));

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j){
            for (int k = 0 ; k < n; ++k){
                AB[i][j] += A[i][k]*B[k][j];
            }
        }
    } 
    return AB;
}

// Norma del suprem matricial
double SupremeNorm(double** A, int n){
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

// Norma sub 1 matricial
double Norm1(double** A, int n){
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

// Escriu una matriu
void PrintMatrix(double** A, int n){
    for (int i = 0 ; i < n; ++i){
        for (int j = 0; j < n; ++j){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl << endl;
    return;
}

// Calcula matriu x vector (retorna vector)
double *MATxVEC(double **M, double x[], int n) {
    double *fx;
    fx = (double*) calloc(n, sizeof(double));
    double num = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            num += M[i][j]*x[j];
        }
        fx[i] = num;
        num = 0;
    }
    return fx;
}

// Norma sub 1 vectorial
double norma1(double v[], int n) {
    double sum = 0; 
    for (int i = 0; i < n; ++i) {
        sum += abs(v[i]);
    }
    return sum;
}
   

// Norma sub infinit vectorial
double norminf(double v[], int n) {
	double maxim = -1;
	for (int i = 0; i < n; ++i) {
		maxim = max(maxim, abs(v[i]));
	}
	return maxim;
}
