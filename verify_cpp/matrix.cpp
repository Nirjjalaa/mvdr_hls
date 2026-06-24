#include "matrix.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

// ================================================
// LOAD 1D CSV
// ================================================
void load_csv_1d(const string& filename, double arr[], int size)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "ERROR: Cannot open " << filename << endl;
        exit(1);
    }

    string line;
    int i = 0;
    while (getline(file, line) && i < size) {
        if (line.empty()) continue;
        // remove any carriage return
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        arr[i++] = stod(line);
    }
    file.close();
}

// ================================================
// LOAD 2D CSV
// ================================================
void load_csv_2d(const string& filename, double arr[][N], int rows, int cols)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "ERROR: Cannot open " << filename << endl;
        exit(1);
    }

    string line;
    int i = 0;
    while (getline(file, line) && i < rows) {
        if (line.empty()) continue;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        stringstream ss(line);
        string val;
        int j = 0;
        while (getline(ss, val, ',') && j < cols) {
            if (!val.empty() && val.back() == '\r')
                val.pop_back();
            arr[i][j++] = stod(val);
        }
        i++;
    }
    file.close();
}

// ================================================
// RECONSTRUCT COMPLEX 1D
// ================================================
void make_complex_1d(double re[], double im[], cplx out[], int size)
{
    for (int i = 0; i < size; i++)
        out[i] = cplx(re[i], im[i]);
}

// ================================================
// RECONSTRUCT COMPLEX 2D
// ================================================
void make_complex_2d(double re[][N], double im[][N], cplx out[][N], int size)
{
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            out[i][j] = cplx(re[i][j], im[i][j]);
}

// ================================================
// VECTOR ERROR
// = ||a - b|| / ||a||
// ================================================
double vec_error(cplx a[], cplx b[], int size)
{
    double diff = 0, norm_a = 0;
    for (int i = 0; i < size; i++) {
        diff   += norm(a[i] - b[i]);
        norm_a += norm(a[i]);
    }
    return sqrt(diff) / sqrt(norm_a);
}

// ================================================
// FROBENIUS ERROR
// = ||A - B||_F / ||A||_F
// ================================================
double frob_error(cplx A[][N], cplx B[][N])
{
    double diff = 0, norm_A = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            diff   += norm(A[i][j] - B[i][j]);
            norm_A += norm(A[i][j]);
        }
    return sqrt(diff) / sqrt(norm_A);
}

// ================================================
// DOT PRODUCT: a^H * b
// ================================================
cplx dot_product(cplx a[], cplx b[], int size)
{
    cplx result = 0;
    for (int i = 0; i < size; i++)
        result += conj(a[i]) * b[i];
    return result;
}
