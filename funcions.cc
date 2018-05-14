#include <iostream>
#include <cmath>
#include "matlib.h"
using namespace std;

// Descomposa A = LU, perm es el vector de permutacions.
int lu(double **a, int n, int perm[], double tol){
    
    int p = 0; 		// Nombre de permutacions.
    
    // Passos del metode de Gauss. 
    for (int k = 0; k < n - 1; ++k) {
        
        
        // Cerquem pivot (pivotatge parcial esglaonat).
        double max = -1; 
        int r = -1;
    	for (int i = k; i < n; ++i) {
    		double s = abs(a[i][k]);
    		for (int j = k + 1; j < n; ++j) {
    			if (s < abs(a[i][j])) s = abs(a[i][j]);	
			}
    		s = abs(a[i][k]/s);
    		if (max < s) {
    			max = s;
    			r = i;
    		}
		}

		if (abs(a[r][k]) < tol) return 0;
		if (k != r) {
			swap(a[k], a[r]);					// Canvi: fila r --- fila k.
	        swap(perm[k], perm[r]);				// Guardem el canvi al vector de permutacions.
	        ++p;	
		}

        
        // Posem zeros a la fila i.
        for (int i = k + 1; i < n; ++i) {
        	double m = a[i][k]/a[k][k];		
        	a[i][k] = m;					
        	
            // Modifiquem tots els altres valors de la fila i
            for (int j = k + 1; j < n; ++j) {
                a[i][j] -= m*a[k][j];
            }
		}
    }
    
	// Retornem 1, -1 en funcio de la paritat del nombre de permutacions de files.
    if (p%2 == 0) return 1;
    return - 1;
    
}

// Resol Ax = b.
void resol(double **a, double x[], double b[], int n, int perm[]){
    
	// Ly = p(permutat)
    double y[n];
    y[0] = b[perm[0]];
    for (int i = 1; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += a[i][j]*y[j];
        }
        y[i] = b[perm[i]] - sum;
    }
    
    // Ux = y
    x[n-1] = y[n-1]/a[n-1][n-1];
    for (int i = n - 2; i >= 0; --i){
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += a[i][j]*x[j];
        }
        x[i] = (y[i] - sum)/a[i][i];
    }
    return;
}

// Calcula la inversa d'A resolent n sistemes per LU.
int inversa(double **a, double **inv_a, int n, double *det_a, double tol){
    
    // Declarem el vector de permutacions.
    int perm[n];
    for (int i = 0 ; i < n; ++i) perm[i] = i;
    
    // Apliquem la descomposicio lu.
    int par = lu(a, n, perm, tol);
    if (par == 0) return 0;
    
    // Calcul del determinant.
    double det = 1;
    for (int i = 0; i < n; ++i) {
    	det *= a[i][i];
	}
	det = det*par;
    *det_a = det;
    
    double e[n];
    for (int i = 0; i < n; ++i) e[i] = 0;
    double x[n];
    
    // Resolem n sistemes lineals canviant el terme independent.
    for (int j = 0; j < n; ++j){
        
        e[j] = 1;
        resol(a, x, e, n, perm);
        for(int i = 0; i < n; ++i) inv_a[i][j] = x[i];
        
        for (int i = 0; i < n; ++i) x[i] = 0;
        e[j] = 0;
    }
    
    return 1;
}
