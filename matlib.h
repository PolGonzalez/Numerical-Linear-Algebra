#ifndef MATLIB_H
#define MATLIB_H
 
double **multiplica_matrius(double** A, double** B,int n,int m,int k);
double **Sum(double** A, double** B, int n);
double **Difference(double** A, double** B, int n);
double SupremeNorm(double** A, int n);
double Norm1(double** A, int n);
void PrintMatrix(double** A, int n);
double **Product(double** A, double** B, int n);
double *MATxVEC(double **M, double x[], int n);
double norma1(double v[], int n);
double norminf(double v[], int n);

#endif
