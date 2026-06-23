#include <iostream>
#include <cmath>
#include "constants.h"
#include "steering_vec.h"
#include "matrix.h"

using namespace std;

int main()
{
    string path = CSV_PATH;

    // ================================================
    // STAGE 1: STEERING VECTOR
    // ================================================
    cout << "===== STAGE 1: STEERING VECTOR =====" << endl;

    cplx a_cpp[N];
    steering_vec(AZ_TARGET, EL_TARGET, D, a_cpp);

    double matlab_real[15] = {
         0.1250, -0.0542, -0.0780,
         0.1218, -0.0277, -0.0978,
         0.1125,  0.0002,  0.1250,
        -0.0542, -0.0780,  0.1218,
        -0.0277, -0.0978,  0.1125
    };
    double matlab_imag[15] = {
         0.0000,  0.1126, -0.0977,
        -0.0279,  0.1219, -0.0778,
        -0.0544,  0.1250,  0.0000,
         0.1126, -0.0977, -0.0279,
         0.1219, -0.0778, -0.0544
    };

    int pass_count = 0;
    for (int i = 0; i < 15; i++) {
        double err = abs(a_cpp[i] - cplx(matlab_real[i], matlab_imag[i]));
        if (err < 1e-3) pass_count++;
    }

    double re1d[N], im1d[N];
    double re2d[N][N], im2d[N][N];

    // compute relative error against MATLAB
    cplx m[15];
    for (int i=0; i<15; i++) m[i] = cplx(matlab_real[i], matlab_imag[i]);
    double rel1 = vec_error(m, a_cpp, 15);

    printf("Elements passed : %d/15\n", pass_count);
    printf("Relative error  : %.6e\n", rel1);
    if (pass_count == 15)
        cout << "STAGE 1: PASS ✅\n" << endl;
    else
        cout << "STAGE 1: FAIL ❌\n" << endl;
// ================================================
// STAGE 2: COVARIANCE MATRIX R
// ================================================
cout << "===== STAGE 2: COVARIANCE MATRIX R =====" << endl;

cplx R_matlab[N][N];
load_csv_2d(path+"R_real.csv", re2d, N, N);
load_csv_2d(path+"R_imag.csv", im2d, N, N);
make_complex_2d(re2d, im2d, R_matlab, N);

// Check 1: Hermitian → R[i][j] = conj(R[j][i])
double herm_err = 0;
for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
        herm_err += abs(R_matlab[i][j] - conj(R_matlab[j][i]));
herm_err /= (N*N);
printf("Hermitian error (avg) : %.6e\n", herm_err);

// Check 2: Diagonal must be real and positive
int diag_pass = 0;
for (int i = 0; i < N; i++)
    if (R_matlab[i][i].real() > 0 && fabs(R_matlab[i][i].imag()) < 1e-6)
        diag_pass++;
printf("Positive real diagonal: %d/%d\n", diag_pass, N);

// Check 3: Trace
double trace = 0;
for (int i = 0; i < N; i++)
    trace += R_matlab[i][i].real();
printf("Trace of R            : %.6f\n", trace);

// Check 4: Print top-left 3x3 (just to visually confirm loading)
cout << "\nTop-left 3x3 of loaded R (real part):" << endl;
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
        printf("%8.4f ", R_matlab[i][j].real());
    cout << endl;
}

if (herm_err < 1e-6 && diag_pass == N)
    cout << "\nSTAGE 2: PASS ✅\n" << endl;
else
    cout << "\nSTAGE 2: FAIL ❌\n" << endl;

    
    return 0;
}
