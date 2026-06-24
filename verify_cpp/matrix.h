#ifndef MATRIX_H
#define MATRIX_H

#include "constants.h"
#include <string>
using namespace std;

// CSV loaders
void load_csv_1d(const string& filename, double arr[], int size);
void load_csv_2d(const string& filename, double arr[][N], int rows, int cols);

// Complex reconstruction
void make_complex_1d(double re[], double im[], cplx out[], int size);
void make_complex_2d(double re[][N], double im[][N], cplx out[][N], int size);

// Error metrics
double vec_error(cplx a[], cplx b[], int size);
double frob_error(cplx A[][N], cplx B[][N]);

// Dot product: result = a^H * b
cplx dot_product(cplx a[], cplx b[], int size);

#endif
