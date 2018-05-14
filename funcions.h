#ifndef FUNCIONS_H
#define FUNCIONS_H
 
int lu(double **a, int n, int perm[], double tol);
void resol(double **a, double x[], double b[], int n, int perm[]);
int inversa(double **a, double **inv_a, int n, double *det_a, double tol);

#endif
